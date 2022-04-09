#!/usr/bin/env python3

from daz_lib import *


#optional:
try:
    import ephem
except:
    print('WARNING: no ephem library, the [optional] dusk/dawn time cannot be calculated')

#for iono correction
import nvector as nv
import iri2016
import pyproj




# get daz iono
################### IONOSPHERE 

def extract_iono_full(esds, framespd):
    # estimating the ionosphere - takes long (several hours)
    # also, storing TECS values (i.e. TEC in slant direction, from IRI2016)
    esds['tecs_A'] = 0.0
    esds['tecs_B'] = 0.0
    esds['daz_iono_grad_mm'] = 0.0
    esds['daz_mm_notide_noiono_grad'] = 0.0
    framespd['Hiono'] = 0.0
    framespd['Hiono_std'] = 0.0
    framespd['Hiono_range'] = 0.0
    framespd['tecs_A'] = 0.0
    framespd['tecs_B'] = 0.0
    for frame in framespd['frame']:
        print(frame)
        resolution = framespd[framespd['frame'] == frame]['azimuth_resolution'].values[0] # in metres
        try:
            #daz_iono_with_F2 = calculate_daz_iono(frame, esds, framespd)
            #daz_iono_grad, hionos, tecs_A_master, tecs_B_master = calculate_daz_iono(frame, esds, framespd, method = 'gomba', out_hionos = True, out_tec_master = True)
            daz_iono_grad, hionos, tecs_A_master, tecs_B_master, tecs_A, tecs_B = calculate_daz_iono(frame, esds, framespd, method = 'gomba', out_hionos = True, out_tec_all = True)
            hiono = np.mean(hionos)
            hiono_std = np.std(hionos)
        except:
            print('some error occurred here')
            continue
        esds.at[esds[esds['frame']==frame].index, 'daz_iono_grad_mm'] = daz_iono_grad*resolution*1000
        esds.at[esds[esds['frame']==frame].index, 'tecs_A'] = tecs_A
        esds.at[esds[esds['frame']==frame].index, 'tecs_B'] = tecs_B
        #esds.at[esds[esds['frame']==frame].index, 'daz_iono_with_F2'] = daz_iono_with_F2
        esds.at[esds[esds['frame']==frame].index, 'daz_mm_notide_noiono_grad'] = esds[esds['frame']==frame]['daz_mm_notide'] - esds[esds['frame']==frame]['daz_iono_grad_mm'] #*resolution*1000
        framespd.at[framespd[framespd['frame']==frame].index, 'Hiono'] = hiono
        framespd.at[framespd[framespd['frame']==frame].index, 'Hiono_std'] = hiono_std
        framespd.at[framespd[framespd['frame']==frame].index, 'Hiono_range'] = max(hionos)-min(hionos)
        framespd.at[framespd[framespd['frame']==frame].index, 'tecs_A'] = tecs_A_master
        framespd.at[framespd[framespd['frame']==frame].index, 'tecs_B'] = tecs_B_master
        #esds.at[esds[esds['frame']==frame].index, 'daz_mm_notide_noiono_F2'] = esds[esds['frame']==frame]['daz_mm_notide'] - esds['daz_iono_with_F2']*resolution*1000
    return esds, framespd


#######################################
# step 3 - get daz iono
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
    pha_iono = 4*np.pi*k/c/f0*dTEC  # i expect faster propagation through plasma. but it may be opposite, not tested, just fast written
    daz_iono = razi*PRF/2/np.pi/f0 * pha_iono
    daz_iono_ifslowed = -daz_iono
    #tecovl = (TECs_B1 - TEC_master_B1)/(fH*fH) - (TECs_B2 - TEC_master_B2)/(fL*fL)
    #daz_iono = -2*PRF*k*f0/c/dfDC * tecovl
    return daz_iono


def calculate_daz_iono(frame, esds, framespd, method = 'gomba', out_hionos = False, out_tec_master = False, out_tec_all = False):
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
        #tecovl = (selected_frame_esds['TECS_B'] - tec_B_master)/(fL*fL) - (selected_frame_esds['TECS_A'] - tec_A_master)/(fH*fH)
        #daz_iono = 2*PRF*k*f0/c/dfDC * tecovl
        # 04/2022 - actually the squares seem not needed, based directly on iono2phase (see article):
        tecovl = (selected_frame_esds['TECS_B'] - tec_B_master)/fL - (selected_frame_esds['TECS_A'] - tec_A_master)/fH
        daz_iono = 2*PRF*k/c/dfDC * tecovl
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
        elif out_tec_all:
            # experiment...
            return daz_iono, hionos, tec_A_master, tec_B_master, tecs_A, tecs_B
        else:
            return daz_iono, hionos
    else:
        return daz_iono



# extra step: get hours from dawn, or dusk respectively:
def get_hours_from_dusk_dawn(framespd):
    sun = ephem.Sun()
    location = ephem.Observer()
    framespd['hours_from_dusk_dawn'] = 0.0
    #
    for i, frameta in framespd.iterrows():
        if frameta['centre_time'] == '0':
            continue
        lon = frameta['center_lon']
        lat = frameta['center_lat']
        master_time = pd.to_datetime(str(frameta['master'])+'T'+frameta['centre_time'])
        location.lat = str(lat)
        location.lon = str(lon)
        location.date = master_time
        try:
            a = location.previous_rising(sun).datetime()
            aa = (master_time - a).total_seconds()/60/60
            b = location.next_rising(sun).datetime()
            bb = (master_time - b).total_seconds()/60/60
            c = location.previous_setting(sun).datetime()
            cc = (master_time - c).total_seconds()/60/60
            d = location.next_setting(sun).datetime()
            dd = (master_time - d).total_seconds()/60/60
        except:
            continue
        hours = 9999
        for x in [aa,bb,cc,dd]:
            if abs(x) < abs(hours):
                hours = x
        #if sunminutes <0:
        #    sunrise = location.previous_transit(sun).datetime()
        #bb = master_time.time()
        #timediff = (bb.hour*60 + bb.minute) - (sunrise.hour*60 + sunrise.minute)
        framespd.at[i, 'hours_from_dusk_dawn'] = hours
    return framespd
