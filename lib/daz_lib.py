#!/usr/bin/env python3

# general imports
import pandas as pd
import numpy as np
from scipy.constants import speed_of_light
from scipy.constants import pi
from scipy import signal
from scipy.stats import linregress
from sklearn.linear_model import HuberRegressor

import glob, os

import urllib
import requests
from lxml import html

# load the csvs
def load_csvs(esdscsv = 'esds.csv', framescsv = 'frames.csv', core_init = False):
    framespd = pd.read_csv(framescsv)
    esds = pd.read_csv(esdscsv)
    if 'Unnamed: 0' in esds.columns:
        esds = esds.drop('Unnamed: 0', axis=1)
    if 'version' in esds.columns:
        esds = esds.drop('version', axis=1)
    if 'Unnamed: 0' in framespd.columns:
        framespd = framespd.drop('Unnamed: 0', axis=1)
    if not 'epochdate' in esds.columns:
        esds['epochdate'] = esds['epoch'].copy(deep=True)
    if 'epoch' in esds.columns:
        esds = esds.drop('epoch', axis=1)
    esds['epochdate'] = esds.apply(lambda x : pd.to_datetime(str(x.epochdate)).date(), axis=1)
    if core_init:
        mindate = esds['epochdate'].min()
        #maxdate = esds['epochdate'].max()
        esds, framespd = df_preprepare_esds(esds, framespd, mindate)
        #esds = esds.reset_index(drop=True)
        framespd = framespd.reset_index(drop=True)
    return esds, framespd


# general functions
def rad2mm_s1(inrad):
    #speed_of_light = 299792458 #m/s
    radar_freq = 5.405e9  #for S1
    wavelength = speed_of_light/radar_freq #meter
    coef_r2m = -wavelength/4/np.pi*1000 #rad -> mm, positive is -LOS
    outmm = inrad*coef_r2m
    return outmm


def m2deg(meters, lat = 0):
    return meters / (111.32 * 1000 * np.cos(np.radians(lat)))


def EN2azi(N, E, heading = -169):
    #ascending: around -13
    #descending: around -169 (or around 13, doesn't matter)
    alpha = np.deg2rad(heading)
    #thanks Chris Rollins!!!!
    return E*np.sin(alpha)+N*np.cos(alpha)





#######################################
# step 2 - get solid Earth tides
################### SOLID EARTH TIDES

def get_tide_in_azimuth(lat, lon, hei, azi, time1):
    '''
    unfinished function to use ETERNA model from 
    https://github.com/hydrogeoscience/pygtide
    '''
    #i can get heights from frame metadata
    startdate = time1
    import pygtide
    pt = pygtide.pygtide()
    
    #duration is in hours
    #samplrate is in 1/s
    
    #here the azimuth should be + for look direction, so 
    pt.predict(lat, lon, hei, startdate, duration = 0, samprate=0, statazimut = azi, tidalcompo = 3)
    data = pt.results()
    return data.iloc[0]


#######################################
# step 2 - merge and furnish esds and framespd sets:
#######################################

def merge_tides(esds, framespd, earthtides):
    #calculate daz_tide from N,E
    esds['daz_tide_mm'] = 0.0
    for frame in framespd['frame']:
        print(frame)
        frameta = framespd[framespd['frame'] == frame]
        tiderows_frame = earthtides[earthtides['frame'] == frame]
        for i,row in esds[esds['frame'] == frame].iterrows():
            tiderow = tiderows_frame[tiderows_frame[' epoch'] == row['epoch']]
            if tiderow.empty:
                print('error - no epoch in tides')
                continue
            E = tiderow[' dEtide'].values[0]
            N = tiderow[' dNtide'].values[0]
            heading = frameta['heading']
            daz_tide_mm = EN2azi(N, E, heading)*1000
            esds.at[i,'daz_tide_mm'] = daz_tide_mm
    return esds


def df_preprepare_esds(esdsin, framespdin, firstdate = '', countlimit = 25):
    #basic fixes
    esds = esdsin.copy(deep=True)
    framespd = framespdin.copy(deep=True)
    esds['daz_mm'] = 0.0
    esds['daz_cc_mm'] = 0.0
    esds['years_since_beginning'] = 0.0
    framespd['count_all'] = 0
    firstdatei = firstdate
    for frame, group in esds.groupby('frame'):
        if not firstdate:
            firstdatei = group['epochdate'].min()
        frameta = framespd[framespd['frame'] == frame]
        if frameta.empty:
            print('Warning, frame {} not found in framespd, using defaults'.format(frame))
            azimuth_resolution = 14.0
        else:
            azimuth_resolution = float(frameta['azimuth_resolution'])
        count = group.epochdate.count()
        if count < countlimit:
            print('small number of {} samples in frame '.format(str(count))+frame+' - removing')
            #esds = esds[esds['frame'] != frame]
            esds = esds.drop(esds.loc[esds['frame']==frame].index)
            framespd = framespd.drop(framespd.loc[framespd['frame']==frame].index)
            continue
        #remove median from daz_total
        medianvalue = group['daz_total_wrt_orbits'].median()
        group['daz_total_wrt_orbits'] = group['daz_total_wrt_orbits'] - medianvalue
        #save the median correction values to framespd
        framespd.at[frameta.index, 'daz_median_shift_mm'] = medianvalue*azimuth_resolution*1000
        framespd.at[frameta.index, 'count_all'] = int(count)
        group['daz_mm'] = group['daz_total_wrt_orbits']*azimuth_resolution*1000
        group['daz_cc_mm'] = group['daz_cc_wrt_orbits']*azimuth_resolution*1000
        group['years_since_beginning'] = group['epochdate'] - firstdatei
        group['years_since_beginning'] = group['years_since_beginning'].apply(lambda x: float(x.days)/365.25)
        #get std, after detrending - but no need to save daz_detrended_mm now....
        group['daz_detrended_mm'] = signal.detrend(group['daz_mm'])
        framespd.at[frameta.index, 'daz_mm_std_all'] = np.std(group['daz_detrended_mm'])
        #update esds
        esds.update(group['daz_total_wrt_orbits'])
        esds.update(group['daz_mm'])
        esds.update(group['daz_cc_mm'])
        esds.update(group['years_since_beginning'])
    for frame in framespd['frame']:
        if frame not in esds['frame'].values:
            framespd = framespd.drop(framespd.loc[framespd['frame']==frame].index)
    #got those in mm, so we can remove the px values now
    esds = esds.drop('daz_total_wrt_orbits', axis=1)
    esds = esds.drop('daz_cc_wrt_orbits', axis=1)
    #we also do not need the RSLC3 information
    if 'esd_master' in esds.columns:
        esds = esds.drop('esd_master', axis=1)
    # some more to add - daz_tide:
    if 'daz_tide_mm' in esds.columns:
        esds['daz_mm_notide'] = esds['daz_mm'] - esds['daz_tide_mm']
    return esds, framespd



#######################################
# step 4 - get plate motion model
################### ITRF2014

def get_ITRF_ENU(lat, lon):
    url = "https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion/model"
    # data to be sent to api
    data = {'name':'modelform',
        'id':'modelform',
        'lat':str(lat),
        'lon':str(lon),
        'model':'itrf2014',
        'format':'ascii'}
    # sending post request and saving response as response object
    r = requests.post(url = url, data = data)
    #outputs are in mm/year, first E, then N
    cont = html.fromstring(r.content)
    [E,N] = cont.text_content().split()[14:16]
    E = float(E)
    N = float(N)
    return E, N


#get ITRF averages within the frames
def df_get_itrf_slopes(framespd):
    framespd['slope_plates_vel_azi_itrf2014'] = 0.0
    framespd['slope_plates_vel_azi_itrf2014_point'] = 0.0
    for ind,frameta in framespd.iterrows():
        frame = frameta['frame']
        itrfs = []
        print(frame)
        clon = frameta['center_lon']
        clat = frameta['center_lat']
        heading = frameta['heading']
        try:
            E, N = get_ITRF_ENU(clat, clon)
        except:
            try:
                E, N = get_ITRF_ENU(clat, clon)
            except:
                print('error')
                continue
        itrf_point = EN2azi(N, E, heading)
        # use a median over 'whole' frame:
        for i in range(round(clon*10-23.4/2),round(clon*10+23.4/2)+1,5):
            lon = i/10
            for j in range(round(clat*10-23.4/2),round(clat*10+23.4/2)+1,5):
                lat = j/10
                try:
                    E, N = get_ITRF_ENU(lat, lon)
                    itrfs.append(EN2azi(N, E, heading))
                except:
                    print('connection error')
        vel_plates_azi = np.mean(itrfs)
        print('difference: '+str(vel_plates_azi - itrf_point))
        #
        framespd.at[ind, 'slope_plates_vel_azi_itrf2014'] = vel_plates_azi
        framespd.at[ind, 'slope_plates_vel_azi_itrf2014_point'] = itrf_point
    return framespd


def df_compare_new_orbits(esds):
    std_diffs = []
    for frame, selected_frame_esds in esds.groupby('frame'):
        neworb = selected_frame_esds[selected_frame_esds['epochdate'] > pd.Timestamp('20200731')]
        oldorb = selected_frame_esds[selected_frame_esds['epochdate'] > pd.Timestamp('20170101')]
        oldorb = oldorb[oldorb['epochdate'] < pd.Timestamp('20200730')]
        if (not neworb.empty) and (not oldorb.empty):
            if not len(oldorb) < len(neworb):
                oldorb = oldorb.tail(len(neworb))
                std_old = oldorb['daz_mm_notide_noiono_grad_OK'].std()
                std_new = neworb['daz_mm_notide_noiono_grad_OK'].std()
                std_diffs.append(std_new - std_old)
    std_diffs = np.array(std_diffs)
    return std_diffs
