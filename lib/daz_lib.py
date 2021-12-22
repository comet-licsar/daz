# general imports

import pandas as pd
import numpy as np
from scipy.constants import speed_of_light
from scipy.constants import pi

from scipy import signal
from scipy.stats import linregress
from sklearn.linear_model import HuberRegressor

from LiCSAR_misc import *
import glob, os

#for iono correction
import nvector as nv
import iri2016
import pyproj
import requests
import urllib

#for visualisation and kml export
import hvplot
import holoviews as hv
from holoviews import opts
hv.extension('bokeh')
import simplekml


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
# step 1 - export ESDs to a txt file:
# see $LiCSAR_procdir/esds/esds.sh
# or just get it from licsinfo's esd table..
# (see LiCSquery - get_daz etc.)

#######################################
# step 2 - add additional information to the esds.txt and frames.txt files
# to be run at JASMIN
# this is to prepare the framespd and esds

def s1_azfm(r, t0, azp):
  """azfr = s1_azfm(r, t0, azp)
  Calculate azimuth FM rate given slant range, reference slant-range delay and the azimuth FM rate polynomial for ScanSAR data
  **Arguments:**
  * r:    slant range (meters)
  * t0:   reference slant range time for the polynomial (center swath delay in s)
  * azp:  polynomial coefficients
  **Output:**
  * the function returns the azimuth FM rate"""
  tsr = 2.0 * r / speed_of_light
  dt = tsr - t0
  azfr = azp[0] + dt * (azp[1] + dt*(azp[2] + dt*(azp[3] + dt*azp[4])))
  return azfr


def get_param_gamma(param, parfile, floatt = True, pos = 0):
    a = grep1line(param,parfile).split()[1+pos]
    if floatt:
        a = float(a)
    return a


def generate_framespd(fname = 'esds2021_frames.txt', outcsv = 'framespd_2021.csv'):
    ### fname is input file containing list of frames to generate the frames csv table
    #in the form of:
    # frame,master,center_lon,center_lat
    a = pd.read_csv(fname)
    a['heading']=0.00
    a['azimuth_resolution']=0.00
    a['avg_incidence_angle']=0.00
    a['centre_range_m']=0.00
    a['centre_time']=''
    a['ka']=0.00
    a['kr']=0.00
    a['dfDC'] = 0.00
    for i,row in a.iterrows():
        frame=row['frame']
        #print(frame)
        tr = int(frame[:3])
        metafile = os.path.join(os.environ['LiCSAR_public'], str(tr), frame, 'metadata', 'metadata.txt')
        if not os.path.exists(metafile):
            print('metadata file does not exist for frame '+frame)
            continue
        primepoch = grep1line('master=',metafile).split('=')[1]
        path_to_slcdir = os.path.join(os.environ['LiCSAR_procdir'], str(tr), frame, 'SLC', primepoch)
        if frame == '174A_05407_121212':
            heading = -10.157417
            azimuth_resolution = 13.968690
            avg_incidence_angle = 39.5118
            centre_range_m = 878941.4133
            centre_time = '14:52:00'
    #        kt = 
        else:
            try:
                heading = float(grep1line('heading',metafile).split('=')[1])
                azimuth_resolution = float(grep1line('azimuth_resolution',metafile).split('=')[1])
                avg_incidence_angle = float(grep1line('avg_incidence_angle',metafile).split('=')[1])
                centre_range_m = float(grep1line('centre_range_m',metafile).split('=')[1])
                centre_time = grep1line('center_time',metafile).split('=')[1]
            except:
                print('some error occurred during frame '+frame)
                azimuth_resolution = 0
                avg_incidence_angle = 0
                centre_range_m = 0
                centre_time = 0
                heading = 0
    #        kt = float(grep1line('kt=',metafile).split('=')[1])
        try:
            dfDC, ka, kr = get_dfDC(path_to_slcdir)
        except:
            print('some error occurred during frame '+frame)
            dfDC = 0
            ka = 0
            kr = 0
        a.at[i,'heading'] = heading
        a.at[i,'azimuth_resolution']  = azimuth_resolution
        a.at[i,'avg_incidence_angle']  = avg_incidence_angle
        a.at[i,'centre_range_m']  = centre_range_m
        a.at[i,'centre_time']  = centre_time
    #    a.at[i,'kt']  = kt
        a.at[i,'dfDC']  = dfDC
        a.at[i,'ka']  = ka
        a.at[i,'kr']  = kr
    a.to_csv(outcsv, float_format='%.4f', index=False)

def get_dfDC(path_to_slcdir, f0=5405000500, burst_interval = 2.758277, returnka = True):
    #f0 = get_param_gamma('radar_frequency', parfile)
    #burst_interval = get_param_gamma('burst_interval', topsparfile)
    parfile = glob.glob(path_to_slcdir+'/????????.slc.par')[0]
    topsparfiles = glob.glob(path_to_slcdir+'/????????.IW?.slc.TOPS_par')
    iwparfiles = glob.glob(path_to_slcdir+'/????????.IW?.slc.par')
    #
    lam = speed_of_light / f0
    dfDC = []
    kas = []
    #krs = []
    #print('This is a proper solution but applied to primary SLC image. originally it is applied by GAMMA on the RSLC...')
    for n in range(len(topsparfiles)):
        topsparfile = topsparfiles[n]
        iwparfile = iwparfiles[n]
        az_steering_rate = get_param_gamma('az_steering_rate', topsparfile) # az_steering_rate is the antenna beam steering rate
        r1 = get_param_gamma('center_range_slc', iwparfile)
        #get the satellite velocity
        #midNstate = int(get_param_gamma('number_of_state_vectors', iwparfile)/2)+1
        # ... actually number of burst info differs... so just using the 1st burst - as anyway we do quite drastic change to dfDC - mean from swaths
        midNstate = 1
        sv = 'state_vector_velocity_' + str(midNstate)
        velc1 = get_param_gamma(sv, iwparfile, pos=0)
        velc2 = get_param_gamma(sv, iwparfile, pos=1)
        velc3 = get_param_gamma(sv, iwparfile, pos=2)
        vsat = np.sqrt(velc1**2 + velc2**2 + velc3**2)
        # now some calculations
        afmrate_srdelay = get_param_gamma('az_fmrate_srdelay_'+ str(midNstate), topsparfile)
        afmrate_poly = []
        afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 0))
        afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 1))
        afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 2))
        try:
            afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 3))
        except:
            afmrate_poly.append(0)
        try:
            afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 4))
        except:
            afmrate_poly.append(0)
        ka = s1_azfm(r1, afmrate_srdelay, afmrate_poly) #unit: Hz/s == 1/s^2
        kr = -2.0 * vsat * az_steering_rate*(pi / 180.0) / lam
        if (kr != 0.0):
            #kt = ka * kr/(kr - ka)
            # but maybe should be kt = (kr*ka)/(ka-kr) # see https://iopscience.iop.org/article/10.1088/1755-1315/57/1/012019/pdf  --- and remotesensing-12-01189-v2, and Fattahi et al...
            # ok, gamma reads kr to be ka... corrected
            kt = kr * ka/(ka - kr)
        else:
            kt = -ka
        #finally calculate dfDC:
        #burst_interval = get_param_gamma('burst_interval', topsparfile)
        kas.append(ka)
        #krs.append(kr)
        dfDC.append(kt*burst_interval) #burst_interval is time within the burst... we can also just calculate.. see Grandin: eq 15: hal.archives-ouvertes.fr/hal-01621519/document
        #ok, that's the thing - burst_interval is actually t(n+1) - t(n) - see remotesensing-12-01189-v2
        #so it should be kt * -burst_interval, that is why GAMMA has the -kt J ... ok, good to realise this
    dfDC = np.mean(dfDC)
    ka = np.mean(kas)
    #kr = np.mean(krs)
    if returnka:
        return dfDC, ka #, kr
    else:
        return dfDC



#######################################
# step 3 - get solid Earth tides
################### SOLID EARTH TIDES



#######################################
# step 4 - merge and furnish esds and framespd sets:
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
    if 'Unnamed: 0' in esds.columns:
        esds = esds.drop('Unnamed: 0', axis=1)
    if 'version' in esds.columns:
        esds = esds.drop('version', axis=1)
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
# step 5 - get daz iono
################### IONOSPHERE 

def get_tecs(glat, glon, altitude, acq_times, returnhei = False):
    '''
    here, the altitude is the satellite altitude (max iono height) to check...
    '''
    altkmrange = [0, altitude, altitude]
    TECs = []
    heis = []
    for acqtime in acq_times:
        iri_acq = iri2016.IRI(acqtime, altkmrange, glat, glon )
        TECs.append(iri_acq.TEC.values[0])
        heis.append(iri_acq.hmF2.values[0])
    if returnhei:
        return TECs, heis
    else:
        return TECs


#get satellite position - ECEF
def aer2ecef(azimuthDeg, elevationDeg, slantRange, obs_lat, obs_long, obs_alt):
    '''
    obs_alt must be in metres
    solution by https://stackoverflow.com/questions/15954978/ecef-from-azimuth-elevation-range-and-observer-lat-lon-alt
    '''
    sitex, sitey, sitez = latlonhei2ecef(obs_lat,obs_long,obs_alt)
    #some needed calculations
    slat = np.sin(np.radians(obs_lat))
    slon = np.sin(np.radians(obs_long))
    clat = np.cos(np.radians(obs_lat))
    clon = np.cos(np.radians(obs_long))
    #
    azRad = np.radians(azimuthDeg)
    elRad = np.radians(elevationDeg)
    #
    # az,el,range to sez convertion
    south  = -slantRange * np.cos(elRad) * np.cos(azRad)
    east   =  slantRange * np.cos(elRad) * np.sin(azRad)
    zenith =  slantRange * np.sin(elRad)
    #
    #
    x = ( slat * clon * south) + (-slon * east) + (clat * clon * zenith) + sitex
    y = ( slat * slon * south) + ( clon * east) + (clat * slon * zenith) + sitey
    z = (-clat *        south) + ( slat * zenith) + sitez
    #
    return x, y, z


def latlonhei2ecef(lat, lon, alt):
    '''
    altitude should be in metres!!!!!
    '''
    transformer = pyproj.Transformer.from_crs(
        {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
        {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
        )
    x, y, z = transformer.transform(lon, lat, alt, radians=False)
    return x, y, z


def ecef2latlonhei(x, y, z):
    transformer = pyproj.Transformer.from_crs(
        {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
        {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
        )
    lon, lat, alt = transformer.transform(x,y,z,radians=False)
    return lat, lon, alt


def get_altitude(lat, lon):
    '''
    uses USGS site to get elevation data. thanks to:
    https://gis.stackexchange.com/questions/338392/getting-elevation-for-multiple-lat-long-coordinates-in-python
    '''
    #
    # USGS Elevation Point Query Service
    #url = r'https://nationalmap.gov/epqs/pqs.php?'
    #
    # opentopodata
    #url = r'https://api.opentopodata.org/v1/eudem25m?locations=51.875127,-3.341298
    url = r'https://api.opentopodata.org/v1/etopo1?locations={0},{1}'.format(lat, lon)
    #
    result = requests.get(url) # + urllib.parse.urlencode(params)))
    elev = result.json()['results'][0]['elevation']
    #print(elev)
    if float(elev) < -100:
        elev = 0
    return elev



def get_abs_iono_corr(frame,esds,framespd):
    selected_frame_esds = esds[esds['frame'] == frame].copy()
    frameta = framespd[framespd['frame']==frame]
    PRF = 486.486
    k = 40.308193 # m^3 / s^2
    f0 = 5.4050005e9
    c = speed_of_light
    razi = frameta.azimuth_resolution.values[0]*1000
    selected_frame_esds['TECS'] = (selected_frame_esds['tecs_A'] + selected_frame_esds['tecs_B'])/2
    master_tecs = ((frameta['tecs_A'] + frameta['tecs_B'])/2).values[0]
    dTEC = selected_frame_esds['TECS'] - master_tecs
    pha_iono = -4*np.pi*k/c/f0*dTEC  # i expect faster propagation through plasma. but it may be opposite
    daz_iono = razi*PRF/2/np.pi/f0 * pha_iono
    daz_iono_ifslowed = -daz_iono
    #tecovl = (TECs_B1 - TEC_master_B1)/(fH*fH) - (TECs_B2 - TEC_master_B2)/(fL*fL)
    #daz_iono = -2*PRF*k*f0/c/dfDC * tecovl
    return daz_iono


def calculate_daz_iono(frame, esds, framespd, method = 'gomba', out_hionos = False, out_tec_master = False):
    '''
    use method either 'gomba' - only gradient, or 'liang' that includes also some extra F2 height correction..
    '''
    selected_frame_esds = esds[esds['frame'] == frame].copy()
    frameta = framespd[framespd['frame']==frame]
    # extract some variables
    heading = frameta['heading'].values[0]
    scene_center_lon = frameta['center_lon'].values[0]
    scene_center_lat = frameta['center_lat'].values[0]
    resolution_of_pixel = frameta['azimuth_resolution'].values[0]
    range_avg = frameta['centre_range_m'].values[0]
    master = frameta['master'].values[0]
    inc_angle_avg = frameta['avg_incidence_angle'].values[0]
    center_time=frameta['centre_time'].values[0]
    dfDC = frameta['dfDC'].values[0]
    #a bit recalculate
    theta = np.radians(inc_angle_avg)
    #sathei = int(range_avg * np.cos(theta)/1000) #in km --- will do this better
    master_time = pd.to_datetime(str(master)+'T'+center_time)
    acq_times = pd.to_datetime(selected_frame_esds.epochdate.astype(str)+'T'+center_time)
    # include master time!
    acq_times[acq_times.index[-1]+1] = master_time
    #
    # 1. get f2 hei inbetween target center point C and nadir of the satellite satg
    wgs84 = nv.FrameE(name='WGS84')
    Pscene_center = wgs84.GeoPoint(latitude=scene_center_lat, longitude=scene_center_lon, degrees=True)
    burst_len = 7100*2.758277 #approx. satellite velocity on the ground 7100 [m/s] * burst_interval [s]
    ###### do the satg_lat, lon
    azimuthDeg = heading-90 #yes, azimuth is w.r.t. N (positive to E)
    elevationDeg = 90-inc_angle_avg
    slantRange = range_avg
    try:
        scene_alt = get_altitude(scene_center_lat, scene_center_lon)
    except:
        scene_alt = 0
    #to get position of the satellite - UNCLEAR about slantRange - is this w.r.t. DEM? (should ask GAMMA) - if not (only in WGS-84), i should use scene_alt=0!
    x, y, z = aer2ecef(azimuthDeg, elevationDeg, slantRange, scene_center_lat, scene_center_lon, scene_alt)
    satg_lat, satg_lon, sat_alt = ecef2latlonhei(x, y, z)
    Psatg = wgs84.GeoPoint(latitude=satg_lat, longitude=satg_lon, degrees=True)
    # get middle point between scene and sat - and get F2 height for it
    path = nv.GeoPath(Pscene_center.to_nvector(), Psatg.to_nvector())
    # get point in the middle
    Pmid_scene_sat = path.interpolate(0.5).to_geo_point()
    # get hionos in that middle point:
    tecs, hionos = get_tecs(Pmid_scene_sat.latitude_deg, Pmid_scene_sat.longitude_deg, 800, acq_times, returnhei = True)
    hiono_master = hionos[-1]
    selected_frame_esds['hiono'] = hionos[:-1]  ###*1000 # convert to metres, avoid last measure, as this is 'master'
    # work in dedicated table
    df = pd.DataFrame(acq_times)
    df['hiono'] = hionos
    #
    ############## now calculate TEC using the SLM knowledge, i.e. different A,B per epoch (!)
    # (note that the last hiono is for the master/reference epoch
    tecs_A = []
    tecs_B = []
    #for hiono in hionos:
    for i,a in df.iterrows():
        hiono = a['hiono']*1000 # m
        epochdate = a['epochdate']
        # first, get IPP - ionosphere pierce point
        # range to IPP can be calculated using:
        range_IPP = hiono/np.sin(theta)
        x, y, z = aer2ecef(azimuthDeg, elevationDeg, range_IPP, scene_center_lat, scene_center_lon, scene_alt)
        ippg_lat, ippg_lon, ipp_alt = ecef2latlonhei(x, y, z)
        Pippg = wgs84.GeoPoint(latitude=ippg_lat, longitude=ippg_lon, degrees=True)
        # then get A', B'
        PsatgA, _azimuth = Psatg.displace(distance=burst_len/2, azimuth=heading-180, method='ellipsoid', degrees=True)
        PsatgB, _azimuth = Psatg.displace(distance=burst_len/2, azimuth= heading, method='ellipsoid', degrees=True)
        # then do intersection ...
        PippAt, _azimuth = Pippg.displace(distance=burst_len, azimuth=heading-180, method='ellipsoid', degrees=True)
        PippBt, _azimuth = Pippg.displace(distance=burst_len, azimuth= heading, method='ellipsoid', degrees=True)
        path_ipp = nv.GeoPath(PippAt, PippBt)
        path_scene_satgA = nv.GeoPath(Pscene_center, PsatgA)
        path_scene_satgB = nv.GeoPath(Pscene_center, PsatgB)
        # these two points are the ones where we should get TEC
        PippA = path_ipp.intersect(path_scene_satgA).to_geo_point()
        PippB = path_ipp.intersect(path_scene_satgB).to_geo_point()
        ######### get TECS for A, B
        TECV_A = get_tecs(PippA.latitude_deg, PippA.longitude_deg, round(sat_alt/1000), [epochdate], False)[0]
        TECV_B = get_tecs(PippB.latitude_deg, PippB.longitude_deg, round(sat_alt/1000), [epochdate], False)[0]
        # get inc angle at IPP - see iono. single layer model function
        earth_radius = 6378160 # m
        sin_thetaiono = earth_radius/(earth_radius+hiono) * np.sin(theta)
        TECS_A = TECV_A/np.sqrt(1-sin_thetaiono**2)
        TECS_B = TECV_B/np.sqrt(1-sin_thetaiono**2)
        tecs_A.append(TECS_A)
        tecs_B.append(TECS_B)
    #
    tec_A_master = tecs_A[-1]
    tec_B_master = tecs_B[-1]
    tecs_A = tecs_A[:-1]
    tecs_B = tecs_B[:-1]
    #
    selected_frame_esds['TECS_A'] = tecs_A
    selected_frame_esds['TECS_B'] = tecs_B
    #
    ##############################
    #if method == 'gomba':
    PRF = 486.486
    k = 40.308193 # m^3 / s^2
    f0 = 5.4050005e9
    c = speed_of_light
    fH = f0 + dfDC*0.5
    fL = f0 - dfDC*0.5
    #tecovl = (TECs_B1 - TEC_master_B1)/(fH*fH) - (TECs_B2 - TEC_master_B2)/(fL*fL)
    #daz_iono = -2*PRF*k*f0/c/dfDC * tecovl
    if method == 'gomba':
        # 08/2021 - empirically checked, correct:
        tecovl = (selected_frame_esds['TECS_B'] - tec_B_master)/(fL*fL) - (selected_frame_esds['TECS_A'] - tec_A_master)/(fH*fH)
        daz_iono = 2*PRF*k*f0/c/dfDC * tecovl
    else:
        # following the Liang 2019:
        # needs ka...
        hsat =  sat_alt
        ka = frameta['ka'].values[0]
        #print('setting Hiono={} km, but it would be close to Gomba method as around 700 km'.format(hiono_m))
        # kion should be rad/s
        # this would get iono effect within the full burst spectrum (from fH to fL)
        phaionL = -4*np.pi*k/c/fL * selected_frame_esds['TECS_B']  # - TEC_master_B1)
        phaionH = -4*np.pi*k/c/fH * selected_frame_esds['TECS_A']  # - TEC_master_B2)
        phaionL_m = -4*np.pi*k/c/fL * tec_B_master
        phaionH_m = -4*np.pi*k/c/fH * tec_A_master
        #
        #phi_ionoramp_burst = fL*fH/(f0*(fH*fH - fL*fL)) * (dphaionL*fH - dphaionH*fL)
        phi_ionoramp_burst = phaionH - phaionL
        phi_ionoramp_burst_m = phaionH_m - phaionL_m
        #i suppose i relate it to burst overlap interval (i.e. 2.7 s /10)
        #burst_ovl_interval = 2.758277/10 #s
        burst_interval = 2.758277
        #kion = phi_ionoramp_burst / burst_ovl_interval # rad/s
        kion = phi_ionoramp_burst / burst_interval # rad/s
        kion_m = phi_ionoramp_burst_m / burst_interval # rad/s
        nesd_iono = -(1/(2*np.pi)) * kion/ka * ((selected_frame_esds['hiono'])/hsat + 0.5) # only for one acq!!!
        nesd_iono_master = -(1/(2*np.pi)) * kion_m/ka * ((hiono_master)/hsat + 0.5) # only for one acq!!!
        # in seconds, need to convert to pixels
        #phiesd = nesd_iono * 2 * np.pi * dfDC
        #daz_iono = PRF * phiesd / (2*np.pi*dfDC)
        daz_iono = PRF * (nesd_iono - nesd_iono_master) # / (2*np.pi*dfDC)
        #return daz_iono
    if out_hionos:
        if out_tec_master:
            return daz_iono, hionos, tec_A_master, tec_B_master
        else:
            return daz_iono, hionos
    else:
        return daz_iono



#######################################
# step 6 - get plate motion model
################### ITRF2014

import requests
from lxml import html

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







#######################################
# step 7 - estimate velocities (and remove outliers etc.)
###################
def get_rmse(y, y_pred, ddof=2):
    return np.sqrt( np.sum(np.square(y - y_pred)) / (len(y)-ddof) )




def get_stdvel(rmse, tmatrix):
    count = tmatrix.shape[0]
    #add column of ones to the tmatrix:
    cones = np.ones(tmatrix.shape)
    A = np.append(tmatrix, cones, axis = 1)
    #
    #now add the rmse to the diagonal of var-covar matrix, already invert it for later
    Qd = np.zeros((count,count)) #,len(d)))
    np.fill_diagonal(Qd,1/rmse**2)
    #
    # do Qm = (G.T Qd^-1 G)^-1  ---- Qm = [var_vel, var_intercept]
    Qm = np.linalg.inv(A.transpose() @ Qd @ A)
    var_vel = Qm.diagonal()[0]
    #var_intercept = Qm.diagonal()[1]
    STD_vel = np.sqrt(var_vel)
    #STD_intercept = np.sqrt(var_intercept)
    return STD_vel


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


def dates2ordinal(x):
    toord = lambda x: x.toordinal()
    vfunc = np.vectorize(toord)
    return vfunc(x)


def fit_huber(x, y, alpha = 1, epsilon=1.2, outinyears = True):
    '''
    just a Huber regression helper function - in case input x is in date(time),
    you may set either the output will be in years (e.g. mm/year) or if not, in days.
    '''
    typex = str(type(x[0])).split("'")[1].split('.')[0]
    xisdate = False
    if typex == 'datetime' or typex == 'pandas':
        xisdate = True
        #it is datetime, so converting to ordinals
        x = dates2ordinal(x)
        x = x - x[0]
    X = np.array([x]).transpose()
    huber = HuberRegressor(alpha = alpha, epsilon = epsilon)
    huber.fit(X, y)
    slope = huber.coef_[0]
    intercept = huber.intercept_
    y_pred = huber.predict(X)
    rmse = get_rmse(y,y_pred)
    std_vel = get_stdvel(rmse, X)
    if xisdate and outinyears:
        slope = slope*365.25
        std_vel = std_vel*365.25
    return slope, intercept, std_vel, ypred, huber.outliers_


def df_calculate_slopes(esdsin, framespdin, alpha = 2.5, eps = 1.5, bycol = 'daz_mm_notide', subset = True, roll_assist = False):
    esds = esdsin.copy(deep=True)
    framespd = framespdin.copy(deep=True)
    #must start with True to initialise all 'as outliers'
    esds['is_outlier_'+bycol] = True
    framespd[bycol+'_RMSE_selection'] = 1.0
    framespd[bycol+'_count_selection'] = 1
    framespd[bycol+'_RMSE_full'] = 1.0
    framespd['slope_'+bycol+'_mmyear'] = -999
    framespd['intercept_'+bycol+'_mmyear'] = -999
    #now calculate per frame
    for frame, group in esds.groupby('frame'):
        print(frame)
        frameta = framespd[framespd['frame'] == frame]
        if frameta.empty:
            print('frame data is empty, skipping')
            continue
        # make preselection
        grsel = group.copy(deep=True)
        # get the full dataset, as we will filter it further on
        #limiting the dataset here, as data after 2020-07-30 have way different PODs (v. 1.7)
        #grsel = grsel[grsel['epochdate'] < pd.Timestamp('2020-07-30')] 
        # 2021-10-12 - finally corrected, probably, using daz_ARP - so using full dataset now!!!
        if subset:
            #limiting the dataset here, as often pre-2016 data are influenced by ionosphere (probably)
            grsel = grsel[grsel['epochdate'] > pd.Timestamp('2016-01-01')]  #[grsel['epoch'] < 20200601]
        if len(grsel) < 20:
            print('most of data removed - cancelling for this frame')
            continue
        correctback = np.median(grsel['daz_mm'])
        grsel['daz_mm'] = grsel['daz_mm'] - correctback #np.median(grsel['daz_mm'])
        grsel = grsel[grsel['daz_mm'] > -400][grsel['daz_mm'] < 400]
        if len(grsel) < 15:
            print('most of data removed - cancelling for this frame')
            continue
        # center dataset by mean (mainly to avoid need of intercept - itrf estimate)
        x = grsel['years_since_beginning'].values
        y = grsel[bycol].values
        #center around zero as origin
        meany = np.mean(y)
        meanx = np.mean(x)
        X = np.array([x]).transpose()
        if roll_assist:
            #
            # prepare rolling (median 5) mean-centered data:
            roll = pd.DataFrame(y, index=x)
            roll = roll.rolling(5).mean()
            y2 = roll[0]
            y2 = y2[y2.notna()]
            x2 = y2.index
            meanx2 = meanx + np.mean(x2)
            meany2 = meany + np.mean(y2)
            X2 = x2.values.reshape(len(x2),1)
        # prepare first iteration to huber (warm start) using itrf slope
        #print(frameta)
        itrfslope = frameta.slope_plates_vel_azi_itrf2014.values[0]
        itrfintercept = np.mean(x)*itrfslope * -1 #+ np.median(y)
        xitrf = np.arange(-5,5)
        yitrf = xitrf * itrfslope + itrfintercept
        Xitrf = xitrf.reshape((10, 1))
        # do huber regression on the original selection
        huber = HuberRegressor(alpha = alpha, epsilon = eps, warm_start=True)
        huber.fit(Xitrf, yitrf)
        huber.fit(X, y)
        slope = huber.coef_[0]
        intercept = huber.intercept_ #??
        y_pred = huber.predict(X)
        y_pred_all = y_pred
        if roll_assist:
            # do additional huber regression on the rolling dataset
            huber_roll = HuberRegressor(alpha = alpha, epsilon = eps, warm_start=True)
            huber_roll.fit(Xitrf, yitrf)
            huber_roll.fit(X2, y2)
            slope_roll = huber_roll.coef_[0]
            #intercept_roll = (meanx2 * slope_roll + huber_roll.intercept_) + meany2
            intercept_roll = huber_roll.intercept_   # CHECK IT
            # to predict, use all dates that go to the pre-rolling dataset...
            y_pred_roll = huber_roll.predict(X)
            # if huber roll is closer to itrf, use this
            if abs(slope-itrfslope) > abs(slope_roll-itrfslope):
                print('using rolling value '+str(slope_roll)+' instead of '+str(slope))
                print('(itrf shows '+str(itrfslope)+' )')
                slope = slope_roll
                intercept = intercept_roll
                y_pred = y_pred_roll
        # calculate RMSE using all points:
        rmse_full = get_rmse(y,y_pred_all)
        std_vel_full = get_stdvel(rmse_full, X)
        # and use outliers from huber !!!!!! even if using huber_rolling !!!!
        grsel['is_outlier_'+bycol] = huber.outliers_
        x = grsel[grsel['is_outlier_'+bycol] == False]['years_since_beginning'].values
        X = np.array([x]).transpose()
        y = grsel[grsel['is_outlier_'+bycol] == False][bycol].values
        try:
            y_pred = huber.predict(X)
            rmse = get_rmse(y,y_pred)
        except:
            print('error getting RMSE (all points are outliers?), using NaN value')
            rmse = np.nan
        #
        framespd.at[frameta.index, 'slope_'+bycol+'_mmyear'] = slope
        framespd.at[frameta.index, 'intercept_'+bycol+'_mmyear'] = intercept + correctback
        # add the rest:
        #keeping the older approach, but perhaps not needed...?
        residuals = y - y_pred
        framespd.at[frameta.index, bycol+'_RMSE_selection'] = rmse
        framespd.at[frameta.index, bycol+'_count_selection'] = len(residuals)
        framespd.at[frameta.index, bycol+'_RMSE_full'] = rmse_full
        # the STD mmy is perhaps the best estimate - kudos to Andy Hooper for this error propagation calculation approach
        framespd.at[frameta.index, bycol+'_RMSE_mmy_full'] = std_vel_full
        #
        # adding std/rmse for velocity - error propagation - this would be ok, but we do not know exact/real shift - so using based on the model here!
        years = grsel.years_since_beginning.max() - grsel.years_since_beginning.min()
        shift_mm = slope * years
        # std vel is from http://ipl.physics.harvard.edu/wp-uploads/2013/03/PS3_Error_Propagation_sp13.pdf - 'should' be right, but perhaps it is not?
        rms_vel_full = abs((rmse_full/shift_mm) * slope)
        rms_vel_sel = abs((rmse/shift_mm) * slope)
        framespd.at[frameta.index, bycol+'_RMSE_mmy_full_error_multi'] = rms_vel_full
        esds.update(grsel['is_outlier_'+bycol])
    return esds, framespd


