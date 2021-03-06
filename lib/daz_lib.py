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

import geopandas
import shapely

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
    lenframes = len(framespd['frame'])
    # oh.. would be better using framespd.iterrows but ok..
    for i, frame in framespd['frame'].iteritems():
        print('  Running for {0:6}/{1:6}th frame...'.format(i+1, lenframes), flush=True, end='\r')
        #print(frame)
        try:
            frameta = framespd[framespd['frame'] == frame]
            tiderows_frame = earthtides[earthtides['frame'] == frame]
            for i,row in esds[esds['frame'] == frame].iterrows():
                tiderow = tiderows_frame[tiderows_frame[' epoch'] == row['epoch']]
                if tiderow.empty:
                    print('error in frame '+frame+'- no epoch in tides')
                    continue
                E = tiderow[' dEtide'].values[0]
                N = tiderow[' dNtide'].values[0]
                heading = frameta['heading']
                daz_tide_mm = EN2azi(N, E, heading)*1000
                esds.at[i,'daz_tide_mm'] = daz_tide_mm
        except:
            print('some error for frame '+frame)
    print('\ndone')
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
################### ITRF2014 (or GSRM_2014)

def get_ITRF_ENU(lat, lon, model='itrf2014', refto='NNR'):
    '''Gets plate motion model values from UNAVCO website
    
    Args
        lat (float): latitude
        lon (float): longitude
        model (string): choose model code, e.g. 'gsrm_2014', 'itrf2014',.. see lines after l. 571 of `view-source:https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion-calculator.html`
        refto (string): choose reference plate, e.g. 'NNR', 'EU',... see lines on the link above after line 683
    
    Returns
        float, float: east and north component of plate motion in given coordinates
    '''
    url = "https://www.unavco.org/software/geodetic-utilities/plate-motion-calculator/plate-motion/model"
    # data to be sent to api
    data = {'name':'modelform',
        'id':'modelform',
        'lat':str(lat),
        'lon':str(lon),
        'model':model,
        'reference':refto,
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


def df_compare_new_orbits(esds, col = 'daz_mm_notide_noiono_grad_OK'):
    '''
    first attempt, not really used function
    '''
    std_diffs = []
    for frame, selected_frame_esds in esds.groupby('frame'):
        neworb = selected_frame_esds[selected_frame_esds['epochdate'] > pd.Timestamp('20200730')]
        oldorb = selected_frame_esds[selected_frame_esds['epochdate'] > pd.Timestamp('20170101')]
        oldorb = oldorb[oldorb['epochdate'] < pd.Timestamp('20200730')]
        if (not neworb.empty) and (not oldorb.empty):
            if not len(oldorb) < len(neworb):
                oldorb = oldorb.tail(len(neworb))
                std_old = oldorb[col].std()
                std_new = neworb[col].std()
                std_diffs.append(std_new - std_old)
    std_diffs = np.array(std_diffs)
    return std_diffs




#function for decomposition (by LS inversion)
def decompose_azi2NE(df, col = 'daz_mm_notide_noiono_grad'):
    velcol = 'slope_'+col+'_mmyear'
    rmscol = col+'_RMSE_mmy_full'
    if not velcol in df.columns:
        #print('workaround for a single column')
        velcol = col
        rmscol = col
    A = []
    d = []
    Qt = []
    for i, row in df.iterrows():
        heading = float(row['heading'])
        At = [np.sin(np.radians(heading)), np.cos(np.radians(heading))]
        Qt.append(1/row[rmscol]**2)
        #b.append(row[rmscol]**2)
        A.append(At)
        d.append(row[velcol])
    A = np.array(A)
    d = np.array(d)
    Q = np.zeros((len(d),len(d)))
    np.fill_diagonal(Q,Qt)
    lstsq = np.linalg.lstsq(A,d, rcond=None)
    # Qm will be variance for V,E and for V,N
    try:
        Qm = np.linalg.inv(A.transpose() @ Q @ A)
    except:
        print('matrix for frames listed below is singular! returning 999999999999 for var')
        print(df['frame'].values)
        Qm = np.zeros([2,2])
        np.fill_diagonal(Qm,999999999)
    Qm = np.abs(Qm.diagonal())
    V_E = lstsq[0][0]
    V_N = lstsq[0][1]
    RMSE_E = np.sqrt(Qm[0])
    RMSE_N = np.sqrt(Qm[1])
    return pd.DataFrame({'V_N': pd.Series(V_N),
                         'V_E': pd.Series(V_E),
                         'RMSE_E': pd.Series(RMSE_E),
                         'RMSE_N': pd.Series(RMSE_N),
                         })


def get_GPS_EN(lat, lon, velnc='velok.nc'):
    '''Gets E,N velocities from external file.
    I prepared the velok.nc file from data available in supplementary files of article:
    Basically, I used the existing velocities and correctly georeferenced them to WGS-84 frame.
    '''
    if not os.path.exists(velnc):
        print('the file with velocities does not exist, exiting')
        return False
    E = float(vels.sel(lon=178.5, lat=85.1, method='nearest').VEL_E.values)
    return E,N


# get ITRF N, E values
def get_itrf_EN(df, samplepoints=3, velnc='velok.nc', refto='NNR'):
    '''Gets EN velocities from ITRF2014 plate motion model (auto-extract from UNAVCO website)
    In case velnc exists, it will be used instead of the model. I know it is bit misleading, as this is not ITRF model. Instead,
    I prepared the velok.nc file from data available in supplementary files of article DOI:10.1002/2014GC005407
    Basically, I used the existing velocities and correctly georeferenced them to WGS-84 frame.
    '''
    usevel = False
    if os.path.exists(velnc):
        print('found velocities nc file - using it instead of ITRF2014')
        vels=xr.open_dataset(velnc)
        usevel = True
    itrfs_N = []
    itrfs_E = []
    itrfs_rms_N = []
    itrfs_rms_E = []
    iii = 0
    fullcount = len(df)
    for ind, row in df.iterrows():
        iii = iii+1
        print('getting ITRF for {0}/{1} cells'.format(iii, fullcount))
        clon = row['centroid_lon']
        clat = row['centroid_lat']
        # use a median over 'whole' frame:
        Es = []
        Ns = []
        for i in range(round(clon*10-23.4/2),round(clon*10+23.4/2)+1,samplepoints):
            lon = i/10
            for j in range(round(clat*10-23.4/2),round(clat*10+23.4/2)+1,samplepoints):
                lat = j/10
                if not usevel:
                    try:
                        E, N = get_ITRF_ENU(lat, lon, refto=refto)
                        Es.append(E)
                        Ns.append(N)
                        #itrfs.append(EN2azi(N, E, heading))
                    except:
                        print('connection error')
                else:
                    E = float(vels.sel(lon=lon, lat=lat, method='nearest').VEL_E.values)
                    N = float(vels.sel(lon=lon, lat=lat, method='nearest').VEL_N.values)
                    Es.append(E)
                    Ns.append(N)
        itrfs_E.append(np.mean(Es))
        itrfs_N.append(np.mean(Ns))
        itrfs_rms_E.append(np.std(Es,ddof=1))
        itrfs_rms_N.append(np.std(Ns,ddof=1))
    df['ITRF_N'] = itrfs_N
    df['ITRF_E'] = itrfs_E
    df['ITRF_RMSE_E'] = itrfs_rms_E
    df['ITRF_RMSE_N'] = itrfs_rms_N
    return df
    



def decompose_framespd(framespd, cell_size = 2.25, crs = "EPSG:4326"):
    '''
    cell_size = 2.25  # this is some ~250x250 km
    '''
    framespd['opass'] = framespd['frame'].str[3]
    gdf = geopandas.GeoDataFrame(framespd, 
                geometry=geopandas.points_from_xy(framespd.center_lon, framespd.center_lat),
                crs=crs)
    # establish a grid
    # projection of the grid
    # total area for the grid
    xmin, ymin, xmax, ymax= gdf.total_bounds
    
    # create the cells in a loop
    grid_cells = []
    centroid_lon = []
    centroid_lat = []
    for x0 in np.arange(xmin, xmax+cell_size, cell_size ):
        for y0 in np.arange(ymin, ymax+cell_size, cell_size):
            # bounds
            x1 = x0-cell_size
            y1 = y0+cell_size
            grid_cells.append( shapely.geometry.box(x0, y0, x1, y1)  )
            centroid_lon.append(x0)
            centroid_lat.append(y0)
    
    grid = geopandas.GeoDataFrame(grid_cells, columns=['geometry'], 
                                     crs=crs)
    grid['centroid_lon'] = centroid_lon
    grid['centroid_lat'] = centroid_lat
    
    # merge framespd and the grid
    merged = geopandas.sjoin(gdf, grid, how='left', op='within')
    
    gridgrouped = merged.groupby('index_right')
    gridagg = gridgrouped.agg(count=('opass', 'count'),
                            opass=('opass', list),
                            centroid_lon=('centroid_lon', 'mean'),
                            centroid_lat=('centroid_lat', 'mean'))
    gridagg = gridagg[gridagg['count'] > 1]
    gridagg = gridagg[gridagg['opass'].str.contains('D', regex=False) & gridagg['opass'].str.contains('A', regex=False)]
    
    #now the gridagg contains only A+D cells
    # 1. reduce merged and grouped:
    for i in merged.index.values:
        if merged.loc[i].index_right not in gridagg.index:
            merged = merged.drop(i)
    
    gridgrouped = merged.groupby('index_right')
    # 2. now do the decomposition
    decomposed = gridgrouped.apply(decompose_azi2NE, 'daz_mm_notide_noiono_grad')
    gridagg['VEL_N_noTI'] = decomposed['V_N'].values
    gridagg['VEL_E_noTI'] = decomposed['V_E'].values
    gridagg['RMSE_VEL_N_noTI'] = decomposed['RMSE_N'].values
    gridagg['RMSE_VEL_E_noTI'] = decomposed['RMSE_E'].values
    
    decomposed = gridgrouped.apply(decompose_azi2NE, 'daz_mm_notide')
    gridagg['VEL_N_noT'] = decomposed['V_N'].values
    gridagg['VEL_E_noT'] = decomposed['V_E'].values
    gridagg['RMSE_VEL_N_noT'] = decomposed['RMSE_N'].values
    gridagg['RMSE_VEL_E_noT'] = decomposed['RMSE_E'].values
    
    '''
    decomposed = gridgrouped.apply(decompose_azi2NE, 'daz_mm')
    gridagg['VEL_N'] = decomposed['V_N'].values
    gridagg['VEL_E'] = decomposed['V_E'].values
    gridagg['RMSE_VEL_N'] = decomposed['RMSE_N'].values
    gridagg['RMSE_VEL_E'] = decomposed['RMSE_E'].values
    '''
    gridagg = gridagg.dropna()
    
    return gridagg



# following functions are to get pod offset - Section 5 of the article:

import pandas as pd
import os

def get_s1b_offsets(esds, framespd, col = 'daz_mm_notide_noiono'):
    offsets = []
    for frame in framespd['frame']:
        #print(frame)
        fpd = framespd[framespd['frame'] == frame]
        epd = esds[esds['frame'] == frame]
        offset = get_s1b_offset(epd, fpd)
        offsets.append(offset)
    framespd['s1ab_offset_mm'] = offsets
    return framespd



from daz_timeseries import get_rmse, model_filter

def get_s1b_offset(epd, fpd, col = 'daz_mm_notide_noiono', fix_pod_offset = True, 
                   split_by_pod = True, fit_offset = False, return_model = False, startfromnoiono = True, mincount = 80 ):
    '''
    epd = selected esd pandas dataframe
    fpd - selected frame pd df
    '''
    if fit_offset and fix_pod_offset:
        print('you do not want to do both..')
        return False
    try:
        epd = epd.set_index(epd.epochdate)
        epd = epd.sort_index()
    except:
        print('')
    if startfromnoiono:
        stdate=pd.Timestamp('2016-07-30').date()
        epd = epd[epd.index>stdate]
        if epd.empty:
            return np.nan
    dazes = epd[col].copy() #.values
    poddate = pd.Timestamp('2020-07-30')
    if fix_pod_offset:
        dazes[dazes.index<poddate]-=39
    epochdates = epd.index.values
    years = epd.years_since_beginning.values
    dazes = dazes.values
    masterdate = pd.Timestamp(fpd.master.values[0])
    mastersat = fpd.s1AorB.values[0]
    isB = flag_s1b(epochdates, masterdate, mastersat)
    if not split_by_pod:
        if fit_offset:
            is_pre = (epochdates<poddate).astype(np.int0)
            A = np.vstack((years,np.ones_like(years),isB, is_pre)).T
            model, stderr = model_filter(A, dazes)
            #model = np.linalg.lstsq(A,dazes, rcond=False)[0]
            cAB = model[2]
            pod_offset = model[3]
            if not return_model:
                return pod_offset
        else:
            # now, we do d = A m, where A is of dt, 1, isB:
            A = np.vstack((years,np.ones_like(years),isB)).T
            model, stderr = model_filter(A, dazes)
            #model = np.linalg.lstsq(A,dazes, rcond=False)[0]
    else:
        # now, we do d = A m, where A is of dt, 1, isB_pre, isB_post:
        isB_pre = (epochdates<poddate).astype(np.int0)*isB
        isB_post = (epochdates>=poddate).astype(np.int0)*isB
        minepochs = 10
        if (np.sum(isB_pre) < minepochs) or (np.sum(isB_post) < minepochs):
            return np.nan
        A = np.vstack((years,np.ones_like(years),isB_pre, isB_post)).T
        model, stderr = model_filter(A, dazes)
        #model = np.linalg.lstsq(A,dazes, rcond=False)[0]
        cAB_pre = model[2]
        cAB_post = model[3]
        cdiff = cAB_post - cAB_pre
        if not return_model:
            return cdiff
    if not return_model:
        #v = model[0]
        #c = model[1]
        c_AB = model[2]
        if c_AB == 0:
            return np.nan
        else:
            return c_AB
    else:
        #get rmse -> stderr, return it
        return model, stderr


def flag_s1b(epochdates, masterdate, mastersat = 'A'):
    if mastersat == 'B':
        masterdate = masterdate + pd.Timedelta('6 days')
    isB = []
    for epoch in epochdates:
        # ok, give +- 1 day tolerance due to midnight issue
        if np.abs(np.mod((epoch - masterdate.date()).days, 12)) <= 1:
            isB.append(0)
        else:
            isB.append(1)
    isB = np.array(isB)
    return isB


def get_pod_offset(dazes, years, thresyears = 4, minsamples = 15):
    r2 = np.ones_like(years)
    r3 = np.ones_like(years)
    r2[years >= thresyears] = 0
    r3[years < thresyears] = 0

    # if too small dataset, cancelling
    if np.sum(r3) < minsamples or np.sum(r2) < minsamples:
        return np.nan
    
    # stack them to get A
    A = np.vstack((years,r2,r3)).T
    model = np.linalg.lstsq(A,dazes, rcond=False)[0]
    offset = model[2] - model[1]
    return offset

'''
# used as:
esds = pd.read_csv('esds_updated_iono_ok_20220324.csv')
esds['epochdate'] = esds.apply(lambda x : pd.to_datetime(str(x.epoch)).date(), axis=1)

firstdatei = pd.Timestamp('2016-07-30')
poddate = pd.Timestamp('2020-07-30')
thresyears = aaa.days / 365.25

esds['years_since_beginning'] = esds['epochdate'] - firstdatei.date()
esds['years_since_beginning'] = esds['years_since_beginning'].apply(lambda x: float(x.days)/365.25)
esds = esds[esds['epochdate'] > firstdatei.date()]
frames = esds['frame'].unique()

offsets = []
for frame in frames:
    dazes = esds[esds['frame'] == frame].daz_mm_notide_noiono_ok.values
    years = esds[esds['frame'] == frame].years_since_beginning.values
    offsets.append(get_pod_offset(dazes, years))

offsets = np.array(offsets)
# filtering - offset would not be larger than... 150 mm?
offsets = offsets[np.abs(offsets)<150]

meanoffset = np.nanmean(offsets)
stderr = np.nanstd(offsets, ddof=3)/np.sqrt(len(offsets))
'''
