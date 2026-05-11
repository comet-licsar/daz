#!/usr/bin/env python3

from daz_lib import *


#optional:
try:
    import ephem
except:
    print('WARNING: no ephem library, the [optional] dusk/dawn time cannot be calculated')

#for iono correction
import nvector as nv
try:
    import iri2020 as iri
except:
    print('iri2020 not found, trying iri2016')
    try:
        import iri2016 as iri
    except:
        print('WARNING: iri2016 not found - please use only CODE corrections and set constant alpha')


import pyproj
import numpy as np
import re
import glob
try:
    import wget
    #import zlib
    from LiCSAR_misc import grep1line
except:
    print('error loading libraries to use of CODE (either wget or LiCSAR_misc). Use of CODE will fail')


# get daz iono
################### IONOSPHERE 

def extract_iono_full(esds, framespd, ionosource = 'iri', use_iri_hei=True):
    """ Full extraction of ionospheric effect from ionosource.
    Note this will create column with the phase advanced effect recalculated to apparent azimuth offset [mm] that has opposite sign.
    Therefore this conforms the GRL article and you can subtract this correction from the original values, as usual.
    
    Args:
        ionosource (str):   either 'iri' or 'code'
        use_iri_hei (bool): estimating F2 peak altitude using IRI (recommended), otherwise the 'valid' height of GIM is used (450 km)
    Returns:
        esds, framespd
    """
    # estimating the ionosphere - takes long (several hours)
    # also, storing TECS values (i.e. TEC in slant direction, from IRI2016)
    esds['tecs_A'] = 0.0
    esds['tecs_B'] = 0.0
    esds['daz_iono_mm'] = 0.0
    #esds['daz_mm_notide_noiono_grad'] = 0.0
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
            if use_iri_hei:
                daz_iono_grad, hionos, tecs_A_master, tecs_B_master, tecs_A, tecs_B = calculate_daz_iono(frame, esds, framespd, method = 'gradient', out_hionos = True, out_tec_all = True, ionosource=ionosource, use_iri_hei=use_iri_hei)
                hiono = np.mean(hionos)
                hiono_std = np.std(hionos)
            else:
                daz_iono_grad, tecs_A_master, tecs_B_master, tecs_A, tecs_B = calculate_daz_iono(frame, esds,
                                                                                                         framespd,
                                                                                                         method='gradient',
                                                                                                         out_hionos=False,
                                                                                                         out_tec_all=True,
                                                                                                         ionosource=ionosource,
                                                                                                         use_iri_hei=use_iri_hei)
                hiono = 450
                hiono_std = 0
        except:
            print('some error occurred extracting TEC(s) here')
            continue
        selesds=esds[esds['frame']==frame].copy()
        # 2023/08: changing sign to keep consistent with the GRL article
        selesds['daz_iono_mm'] = daz_iono_grad*resolution*1000
        selesds['tecs_A'] = tecs_A
        selesds['tecs_B'] = tecs_B
        # skipping the correction here, since daz_mm_notide might not exist/not needed:
        #selesds['daz_mm_notide_noiono_grad'] = selesds['daz_mm_notide'] + selesds['daz_iono_grad_mm'] #*resolution*1000
        esds.update(selesds)
        framespd.at[framespd[framespd['frame'] == frame].index[0], 'Hiono'] = hiono
        if use_iri_hei:
            framespd.at[framespd[framespd['frame']==frame].index[0], 'Hiono_std'] = hiono_std
            framespd.at[framespd[framespd['frame']==frame].index[0], 'Hiono_range'] = max(hionos)-min(hionos)
        framespd.at[framespd[framespd['frame']==frame].index[0], 'tecs_A'] = tecs_A_master
        framespd.at[framespd[framespd['frame']==frame].index[0], 'tecs_B'] = tecs_B_master
        #esds.at[esds[esds['frame']==frame].index, 'daz_mm_notide_noiono_F2'] = esds[esds['frame']==frame]['daz_mm_notide'] - esds['daz_iono_with_F2']*resolution*1000
    return esds, framespd


#######################################
# step 3 - get daz iono
################### IONOSPHERE 

def get_tecs(glat, glon, altitude, acq_times, returnhei = False, source='jpl', alpha = 0.85, returnalpha = False, tecxr = None):
    '''Gets estimated TEC over given point, up to given altitude
    
    Args:
        glat, glon, altitude: coordinates and max 'iono height' to get the TEC values for. Altitude is in km
        acq_times (list of dt.datetime): time stamps to get the TEC over the given point
        returnhei (boolean):  if True, it would return TEC values but also estimated F2 peak heights (from IRI)
        source (str): source of TEC - either 'iri' for IRI2016 model (must be installed), 'code' to autodownload from CODE GIM, 'jpl' to use JPL HRES GIM (would fallback to CODE), or 'jpliri' to drive the JPL GIM using IRI2020
        alpha (float): for CODE only, estimate of ratio of TEC towards 'to the satellite only'. If 'auto', it will estimate it using iri. 0.85 is good value
        returnalpha (bool): if True, will also return alpha param (useful with alpha = 'auto')
        tecxr (None or xr.Dataset): only for GIM source. if provided (as get_vtec_from_code), it would use this rather than re-generate (to save some time)
    '''
    if source == 'jpl':
        source = 'code'  # TODO.. for now, both CODE and JPL will try JPL first..
    if glon > 180:
        glon = glon - 180  # fix round-earth
    if returnhei and source != 'iri':
        print('WARNING, height is estimated only through IRI model, now setting to it')
        source = 'iri'
    if alpha == 'auto':
        getalpha=True
    else:
        getalpha=False
    altkmrange = [0, altitude, altitude]
    TECs = []
    heis = []
    alphas = []
    if source == 'jpliri':
        ftable = get_f107_table()
        print('this was an experiment but seems not worth further works')
    for acqtime in acq_times:
        if source == 'iri':
            iri_acq = iri.IRI(acqtime, altkmrange, glat, glon )
            TECs.append(iri_acq.TEC.values[0])
            heis.append(iri_acq.hmF2.values[0])
            if getalpha:
                iri_acq_gps = iri.IRI(acqtime, [0, 20000, 20000], glat, glon)
                alpha = float(iri_acq.TEC / iri_acq_gps.TEC)
                # print('using alpha of ' + str(alpha))
            alphas.append(alpha)
        elif source == 'code':
            if getalpha:
                iri_acq_gps = iri.IRI(acqtime, [0, 20000, 20000], glat, glon )
                iri_acq = iri.IRI(acqtime, altkmrange, glat, glon )
                alpha = float(iri_acq.TEC/iri_acq_gps.TEC)
                # print('using alpha of '+str(alpha))
            alphas.append(alpha)
            try:
                if type(tecxr) != type(None):
                    tec = get_vtec_from_tecxr(tecxr, acqtime, glat, glon, method='linear')
                else:
                    tec = get_vtec_from_code(acqtime, glat, glon)
            except:
                # CODE data is not available earlier than with 6 months delay.. or more?
                tec = np.nan
                print('No JPL/CODE data for date '+str(acqtime.date())+'. Setting NaN.')
            # decrease the value by some alpha... we expect alpha % of TEC being below the satellite.. should be improved
            #alpha = 0.85
            tec = alpha*tec
            TECs.append(tec)
        elif source == 'jpliri':
            # testing - 2025-10-23
            # import iri2020 as iri
            vtec = get_vtec_from_code(acqtime, glat, glon)
            f107 = get_f107_dt(acqtime, ftable)
            drive_iri_with_gim(glat, glon, acqtime, vtec, f107)
            alphas.append(alpha)
    if returnhei:
        if returnalpha:
            return TECs, heis, alphas
        else:
            return TECs, heis
    else:
        if returnalpha:
            return TECs, alphas
        else:
            return TECs



def download_code_data(acqtime, storedir = '/gws/ssde/j25a/nceo_geohazards/vol1/code_iono'):
    """Downloads Ionospheric TEC data from JPL or CODE."""
    ffound = False
    if not ffound:
        filename = 'jpld'+ acqtime.strftime('%j') + '0.' + acqtime.strftime('%y')+'i.nc.gz' # TODO: check YMD
        url = 'https://sideshow.jpl.nasa.gov/pub/iono_daily/gim_for_research/jpld/' + str(acqtime.year) + '/' + filename
        fullpath = os.path.join(storedir, filename)
        ionix = fullpath[:-3]
        if not os.path.exists(ionix):
            if not os.path.exists(fullpath):
                # download this
                try:
                    wget.download(url, out=storedir)
                    ffound = True 
                except:
                    #print('jpl-hr not exist')
                    ffound = False
            if os.path.exists(fullpath):
                ffound = True
        else:
            ffound = True
    ###if JPL HR-GIM is not available, then try to download CODE data. After, 2024-08-01 we need CODE data. 
    if not ffound:  
        instrings = ['CODG', 'CGIM']
        #lastrings = ['I', 'N']
        for instr in instrings:
            #for lastr in lastrings:
            filename = instr + acqtime.strftime('%j') + '0.' + acqtime.strftime('%y') + 'I.Z'
            url = 'http://ftp.aiub.unibe.ch/CODE/' + acqtime.strftime('%Y') + '/' + filename
            fullpath = os.path.join(storedir, filename)
            ionix = fullpath[:-2]
            if not os.path.exists(ionix):
                if not os.path.exists(fullpath):
                    # download this
                    try:
                        wget.download(url, out=storedir)
                        ffound = True 
                    except:
                        #print('code-oldname not exists')
                        ffound = False
                if os.path.exists(fullpath):
                    ffound = True
            else:
                ffound = True
            if ffound:
                break  #Exit loop if successful download
    # since 12/2022 they changed naming convention to e.g. COD0OPSFIN_20230510000_01D_01H_GIM.INX.gz
    # see https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html#iono
    if not ffound:
        filename = 'COD0OPSFIN_'+ acqtime.strftime('%Y') + acqtime.strftime('%j')+'0000_01D_01H_GIM.INX.gz' # TODO: check YMD
        #filename = instr + acqtime.strftime('%j') + '0.' + acqtime.strftime('%y') + 'I.Z'
        url = 'http://ftp.aiub.unibe.ch/CODE/' + acqtime.strftime('%Y') + '/' + filename
        fullpath = os.path.join(storedir, filename)
        ionix = fullpath[:-3]
        if not os.path.exists(ionix):
            if not os.path.exists(fullpath):
                # download this
                try:
                    wgotfile = wget.download(url, out=storedir)
                except:
                    print('error during wget download')
                    ffound = False
            if os.path.exists(fullpath):
                ffound = True
        else:
            ffound = True
    if not ffound:
        print('no CODE layer found for '+filename)
        return False
    if not os.path.exists(ionix):
        rc = os.system('cd ' + storedir + '; 7za x ' + filename + ' >/dev/null 2>/dev/null; rm ' + fullpath)
    if not os.path.exists(ionix):
        print('ERROR: maybe you do not have 7za installed')
        return False
    return ionix


def get_vtec_from_code(acqtime, lat = 0, lon = 0,
                       storedir = '/gws/ssde/j25a/nceo_geohazards/vol1/code_iono', return_fullxr = False, noJPL=False, printout=True):
    """ Adapted from Reza Bordbari script, plus using functions from https://notebook.community/daniestevez/jupyter_notebooks/IONEX
    
    17/03/2025-(MN):function also helps to extract NASA JPL High Resolution vTEC values (15min, 1x1degree) at first. 
    https://sideshow.jpl.nasa.gov/pub/iono_daily/gim_for_research/jpli/
    
    Args:
        acqtime (dt.datetime)
        lat (float)
        lon (float)
        storedir (str)
        return_fullxr (bool): if True, will return full TEC datacube
    """
    #D = acqtime.strftime('%Y%m%d')
    #ipp = np.array([lat,lon])
    # check if exists:

    if not noJPL:
        fna = glob.glob(storedir + '/jpld' + acqtime.strftime('%j') + '0.' + acqtime.strftime('%y') + '*.nc')  # prioritize JPL-HR GIM
        if fna:  
            ionix = os.path.join(storedir, fna[0])  # Found JPL-HR GIM file, use it
        else:
            # JPL-HR GIM does not exist, try to download it.
            ionix = download_code_data(acqtime, storedir)
            fna = glob.glob(storedir + '/jpld' + acqtime.strftime('%j') + '0.' + acqtime.strftime('%y') + '*.nc')  # prioritize JPL-HR GIM again
            if fna:  
                ionix = os.path.join(storedir, fna[0])  # Found a different GIM file, use it
            else:
                ionix = None
    else:
        ionix = None

    if not ionix:
        # If JPL-HR GIM is missing or noJPL is True, fallback to CODE GIM
        fna = glob.glob(storedir + '/????' + acqtime.strftime('%j') + '0.' + acqtime.strftime('%y') + '?')  # CODE GIM
        if fna:
            ionix = os.path.join(storedir, fna[0])
        else:
            # If no CODE GIM is found, try to download it
            ionix = download_code_data(acqtime, storedir)
        if not ionix:
            if printout:
                print('no GIM data available for '+str(acqtime))
            return False
        elif printout:
            print('using CODE GIM data')
    elif printout:
        print('using JPL-HR GIM data')
    #
    #else:
    #    rc=os.system('rm '+fullpath) # clean the .Z
    # prep 
    #hhmmss=acqtime.strftime('%H%M%S')
    # loading the TEC maps, thanks to https://notebook.community/daniestevez/jupyter_notebooks/IONEX (but improved towards xarray by ML B-)
    if os.path.basename(ionix).startswith('jpl'):
        # Open the NetCDF file
        ds = xr.open_dataset(ionix)
        # Convert time epochs to readable datetime format, 15min resolution referenced to j2000 (1/1/2000 12:00 UT).
        time_values = ds['time'].values
        #converting yyyy-mm-ddThh:mm:ss format. #TODO: We can change that regarding ML's wishes. I found this is clear for now.
        converted_times = pd.to_datetime(time_values, origin='2000-01-01 12:00:00', unit='s')
        #dataset2dataarray
        tecxr_jhr = xr.DataArray(data=ds['tecmap'].values, dims=['time','lat','lon'],
                            coords=dict(time=converted_times, lat= ds["lat"].values, lon= ds["lon"].values) )
        tecxr=tecxr_jhr*1e+16 # from TECU   
    else:
        try:
            tecmaps = get_tecmaps(ionix)
        except:
            print('ERROR loading ionix file: '+ionix)
            return False
        try:
            interval = int(grep1line('INTERVAL',ionix).split()[0])
        except:
            print('ERROR, the ionix file '+ionix+' does not contain necessary keywords. Cancelling')
            return False
        timestep = interval/3600
        timecoords = np.arange(0.0,24.0+timestep,timestep)  # we expect start/end time being midnight, should be standard for all CODE files?
        timecoords = acqtime.normalize() + pd.to_timedelta(timecoords, unit='h')  # set to real datetime..
        lat_all = np.arange(87.5,-87.5-2.5,-2.5)
        lon_all = np.arange(-180.0,180.0+5,5.0)
        tecxr = xr.DataArray(data=tecmaps, dims=['time','lat','lon'],
                            coords=dict(time=timecoords, lon=lon_all, lat=lat_all) )
        # interpolate through the nan values
        tonan=9999
        tecxr.where(tecxr!=tonan)
        tecxr=tecxr.interpolate_na(dim="lon", method="linear", fill_value="extrapolate")
        tecxr = tecxr*1e+16 # from TECU
    if acqtime<tecxr.time.min():
        # need to join day before:
        tecxr2=get_vtec_from_code(acqtime-pd.Timedelta('4 hours'), lat=lat, lon=lon, return_fullxr = True)
        tecxr = xr.concat([tecxr2, tecxr], dim='time')
    if acqtime>tecxr.time.max():
        tecxr2=get_vtec_from_code(acqtime+pd.Timedelta('4 hours'), lat=lat, lon=lon, return_fullxr = True)
        tecxr = xr.concat([tecxr, tecxr2], dim='time')
    if return_fullxr:
        return tecxr
    else:
        return get_vtec_from_tecxr(tecxr, acqtime, lat, lon)


# get_vtec_from_code(acqtime, lat, lon, storedir = '/gws/ssde/j25a/nceo_geohazards/vol1/code_iono', return_fullxr = False):
def get_vtec_from_tecxr(tecxr, acqtime, lat, lon, rotate=True, method='cubic'):
    '''
    This will get VTEC from the xr TEC dataset for given acqtime (dt.datetime or pd.Timestamp), lat, lon.
    It is recommended to use 'rotate' that does the 3D interpolation with rotation based on:
      https://github.com/insarlab/MintPy/blob/main/src/mintpy/objects/ionex.py
    that is actually based on
    Schaer, S., Gurtner, W., & Feltens, J. (1998). IONEX: The ionosphere map exchange format
            version 1.1. Paper presented at the Proceedings of the IGS AC workshop, Darmstadt, Germany.
    '''
    # if len(tecxr.time.values) == 25:
    #     # print('old CODE format, 25 values')
    #     h_time = float(acqtime.strftime('%H'))
    #     m_time = float(acqtime.strftime('%M'))
    #     s_time = float(acqtime.strftime('%S'))
    #     # given time in decimal format
    #     time_dec = h_time + (m_time/60) + (s_time / 3600)
    #     # ML: 2023/08, based on : https://github.com/insarlab/MintPy/blob/main/src/mintpy/objects/ionex.py
    #     # that is actually based on
    #     # Schaer, S., Gurtner, W., & Feltens, J. (1998). IONEX: The ionosphere map exchange format
    #     #         version 1.1. Paper presented at the Proceedings of the IGS AC workshop, Darmstadt, Germany.
    #     if rotate:
    #         # 3D interpolation with rotation as above reference
    #         try:
    #             htimes = tecxr.time.values
    #         except:
    #             htimes = np.array([float(tecxr.time)])  #in case of only one value
    #         pretime = int(htimes[htimes <= time_dec][-1])
    #         postime = int(htimes[htimes >= time_dec][0])
    #         #
    #         lon0 = lon + (time_dec - pretime) * 360. / 24.
    #         lon1 = lon + (time_dec - postime) * 360. / 24.
    #         # print(time_dec - pretime)
    #         # print(time_dec - postime)
    #         #
    #         tec_val0 = float(tecxr.interp(time=pretime, lon=lon0, lat=lat, method=method))
    #         tec_val1 = float(tecxr.interp(time=postime, lon=lon1, lat=lat, method=method))
    #         #
    #         tec = ((postime - time_dec) / (postime - pretime) * tec_val0
    #                    + (time_dec - pretime) / (postime - pretime) * tec_val1)
    #     else:
    #         # previous attempt, but still too different from the S1_ETAD CODE outputs (that rotates the Earth towards the Sun..)
    #         tec = float(tecxr.interp(time=time_dec, lon=lon,lat=lat, method=method)) # should be better than linear, but maybe quadratic is more suitable?
    # elif len(tecxr.time.values) > 90: ##Normally tecxr should include each 15min data but, data for the following days cover from 00:00UT to 23:30UT only:6 - 9 Jan 2023, 11 Jan 2023, 13 - 15 Jan 2023, 17 - 18 Jan 2023, 22 Jan 2023.
    #     # print('JP-HR GIM format, 96 values')
    if rotate:
        # 3D interpolation with rotation as above reference
        try:
            htimes = tecxr.time.values
        except:
            htimes = np.array([float(tecxr.time)])  #in case of only one value
        # print(htimes)
        pretime = htimes[htimes <= acqtime].max()
        postime = htimes[htimes >= acqtime].min()
        #convet to timestamps to play with total_seconds
        pretime = pd.Timestamp(pretime)
        postime = pd.Timestamp(postime)

        lon0 = lon + (acqtime - pretime).total_seconds() / 86400 * 360. #let's play with seconds 24x60x60(seconds per day)
        lon1 = lon + (acqtime - postime).total_seconds() / 86400 * 360. #TODO or postime-acqtime?

        # #
        tec_val0 = float(tecxr.interp(time=pretime, lon=lon0, lat=lat, method=method))
        tec_val1 = float(tecxr.interp(time=postime, lon=lon1, lat=lat, method=method))

        tec = ((postime - acqtime).total_seconds() / (postime - pretime).total_seconds() * tec_val0
               + (acqtime - pretime).total_seconds() / (postime - pretime).total_seconds() * tec_val1)     ##linear in time
    else:
        # If rotation is NOT enabled, perform standard cubic interpolation
        tec = float(tecxr.interp(time=acqtime, lon=lon, lat=lat, method=method))
    return tec


def parse_map(tecmap, exponent = -1):
    tecmap = re.split('.*END OF TEC MAP', tecmap)[0]
    return np.stack([np.fromstring(l, sep=' ') for l in re.split('.*LAT/LON1/LON2/DLON/H\\n',tecmap)[1:]])*10**exponent


def get_tecmaps(filename):
    try:
        exponent = int(grep1line('EXPONENT',filename).split()[0]) # this is exponent of the data
    except:
        print('WARNING, exponent not found in '+filename+'. Perhaps the file is corrupted?')
        exponent = -1
    with open(filename) as f:
        ionex = f.read()
        return [parse_map(t, exponent) for t in ionex.split('START OF TEC MAP')[1:]]


def get_tec(tecmap, lat, lon):
    i = round((87.5 - lat)*(tecmap.shape[0]-1)/(2*87.5))
    j = round((180 + lon)*(tecmap.shape[1]-1)/360)
    return tecmap[i,j]


# not used :
#def julian_to_datetime(julian_date):
#    # Julian date of the Unix epoch (1970-01-01)
#    julian_epoch = 2440587.5
#    # Convert to seconds since Unix epoch
#    seconds = (julian_date - julian_epoch) * 86400.0
#    return datetime.utcfromtimestamp(seconds)

def get_f107_table(url='https://spaceweather.gc.ca/solar_flux_data/daily_flux_values/fluxtable.txt'):
    ftable = pd.read_csv(url, sep=r'\s+', skiprows=[1])
    #ftable['fluxdt'] = ftable['fluxjulian'].apply(julian_to_datetime)
    ftable['fluxdatetime'] = pd.to_datetime(
        ftable['fluxdate'].astype(str) + ftable['fluxtime'].astype(str),
        format='%Y%m%d%H%M%S'
    )
    ftable = ftable.set_index('fluxdatetime')
    ftable = ftable.sort_index()
    return ftable


def get_f107_dt(dtime, ftable=None, colname='fluxadjflux'):
    if type(ftable) == type(None):
        ftable = get_f107_table()
    dtimes_numeric = ftable.index.values.astype('datetime64[ns]').astype(np.int64)
    #np.array([dt.value for dt in ftable.index.values], dtype=np.float64)  # .value gives nanoseconds since epoch
    target_numeric = pd.Timestamp(dtime).value
    # Perform linear interpolation
    return np.interp(target_numeric, dtimes_numeric, ftable[colname].values)

'''
unfinished attempt (tried to get chatgpt-supported answer...)
def drive_iri_with_gim(lat, lon, acqtime, GIM_VTEC, f107):
    """
    Adjusts IRI NmF2 and hmF2 to match GIM VTEC.

    Args:
        lat: Latitude (degrees).
        lon: Longitude (degrees).
        acqtime: pd.timestamp or datetime object.
        GIM_VTEC: GIM VTEC value (TECU).
        f107: F10.7 solar radio flux index.

    Returns:
        IRI electron density profile, adjusted NmF2, adjusted hmF2.
    """
import iri2020 # import runiri
try:
    time = acqtime.to_pydatetime()
except:
    time = acqtime
altkmrange = [0,20000,10]
iri_acq = iri2020.IRI(acqtime, altkmrange, glat, glon )
hmf2 =  float(iri_acq['hmF2'])
maxtec_alt = float(iri_acq.alt_km[iri_acq.ne.argmax()])

timedatetime.datetime(2023, 3, 21, 12)  # UTC datetime
glat = 45.0   # Latitude in degrees
glon = -75.0  # Longitude in degrees
alt_km = 300  # Altitude in km
f107 = 150.0  # Custom F10.7 value
f107a = 150.0 # 81-day average F10.7
ap = 4        # Ap index

iri = runiri(
    dtime=time,
    glat=glat,
    glon=glon,
    alts_km=[alt_km],
    f107=f107,
    f107a=f107a,  # optional, see below
    jf=[1]*50
)

iri = iri2016.IRI # sm.SpacepyModel('iri2016')  # Or whichever version you use
iri['year'] = time.year
iri['month'] = time.month
iri['day'] = time.day
iri['hour'] = time.hour
iri['minute'] = time.minute
iri['sec'] = time.second
iri['glat'] = lat
iri['glon'] = lon
iri['altkm'] = np.arange(80, 1000, 10)  # Altitude range (km)
iri['f107'] = f107
iri['f107a'] = f107 # smoothed F10.7
iri.run()

IRI_TEC_initial = iri_acq['TEC'].item() #extract total TEC from 0 to 2000km altitude (or a high altitude)

NmF2_initial = iri_acq['NmF2'].item()
hmF2_initial = iri_acq['hmF2'].item()

NmF2 = NmF2_initial
hmF2 = hmF2_initial
tolerance = 0.1  # TECU
max_iterations = 10

for i in range(max_iterations):
    NmF2_new = NmF2 * (GIM_VTEC / IRI_TEC_initial)
    NmF2 = NmF2_new
    iri_acq['NmF2'] = NmF2  # Set NmF2 manually

    iri.run() # rerun with new Nmf2

    IRI_TEC = iri['tec'].item()
    difference = GIM_VTEC - IRI_TEC
    print(f"Iteration {i+1}: IRI TEC = {IRI_TEC:.2f}, GIM TEC = {GIM_VTEC:.2f}, Difference = {difference:.2f}")

    if abs(difference) < tolerance:
        print("Convergence achieved!")
        break

    return iri['ne'], NmF2, hmF2


# Example usage:
import datetime
lat = 34.0  # Latitude
lon = -118.0 # Longitude
time = datetime.datetime(2023, 10, 27, 12, 0, 0)  # Example time
GIM_VTEC = 30.0  # Example GIM VTEC value (TECU)
f107 = 150.0  # Example F10.7 index

ne_profile, adjusted_NmF2, adjusted_hmF2 = drive_iri_with_gim(lat, lon, time, GIM_VTEC, f107)

print(f"Adjusted NmF2: {adjusted_NmF2:.2f}")
print(f"Adjusted hmF2: {adjusted_hmF2:.2f}")

#Now you can extract the Ne (electron density profile) from iri.
'''



'''

Reza's code (after attempt to convert mat2python - sorry for ugliness here, ML)
def get_vtec_from_code(yyyymmdd = None,day_of_year = None,hhmmss = None,ipp_lat = None,ipp_lon = None): 
    ipp = np.array([ipp_lat,ipp_lon])
    D = num2str(yyyymmdd)
    fullDate = np.array([[str2double(D(np.arange(1,4+1)))],[str2double(D(np.arange(5,6+1)))],[str2double(D(np.arange(7,8+1)))]])
    url = 'http://ftp.aiub.unibe.ch/CODE/' + string(fullDate(1)) + '/' + 'CODG' + string(day_of_year * 10) + '.' + D(np.arange(3,4+1)) + 'I.Z'
    code_zip = 'CODG' + string(day_of_year * 10) + '.' + D(np.arange(3,4+1)) + 'I.Z'
    ionix = 'CODG' + string(day_of_year * 10) + '.' + D(np.arange(3,4+1)) + 'I'
    try:
        zipFile = websave(code_zip,url)
        # the URL was found if you get here.  If not, it goes to the catch.
        print('SUCCESS: Downloaded %s to %s \n' % (url,zipFile))
        command = sprintf('uncompress %s',code_zip)
        system(command)
    finally:
        pass
    
    # Data acquisition time
    T = num2str(hhmmss)
    h_time = str2double(T(np.arange(1,2+1)))
    m_time = str2double(T(np.arange(3,4+1)))
    s_time = str2double(T(np.arange(5,6+1)))
    # given time in decimal format
    min_temp = m_time / 60
    time_dec = h_time + min_temp + (s_time / 3600)
    time_sec = (h_time * 3600) + (m_time * 60) + s_time
    # all the grid points for time
    time_all = np.arange(0,24+1,1)
    # find the 2 closest time indices
    time_sort = np.array([h_time,h_time + 1])
    # longitude modification
    we = 360 / (24 * 60 * 60)
    
    lon_ipp_1 = ipp(2) + (time_sec - (time_sort(1) * 3600)) * we
    lon_ipp_2 = ipp(2) + (time_sec - (time_sort(2) * 3600)) * we
    ## (3) Time and coordinate interpolation
# (3-1) Find indices of the two adjacent TEC maps in latitude
# find coordinates (lat,lon) of the all grid points
    lat_all = np.arange(87.5,- 87.5+- 2.5,- 2.5)
    # find the 2 closest latitude indices
    __,ind_lat = __builtint__.sorted(np.abs(ipp(1) - lat_all),'ascend')
    ind_lat_int = __builtint__.sorted(ind_lat(np.arange(1,2+1)),'ascend')
    # (3-2) Reading IONEX file
    lat_sort,iLat = __builtint__.sorted(lat_all(ind_lat_int),'descend')
    # title of the Global IONEX file (ex. http://ftp.aiub.unibe.ch/CODE/2020/)
    IONEXFile = ionix
    latBlock1_t1,latBlock2_t1,latBlock1_t2,latBlock2_t2 = get_lat_block(IONEXFile,fullDate,time_dec,time_sort,lat_sort)
    # (3-3) Bilinear interpolation (interpolation in lat, long, time)
# theory is based on: https://www.omnicalculator.com/math/bilinear-interpolation
    
    lon_all = np.arange(- 180,180+5,5)
    # epoch one
# find the 2 closest longitude indices (based on the 1st modified IPP longitude)
    __,ind_lon = __builtint__.sorted(np.abs(lon_ipp_1 - lon_all),'ascend')
    ind_lon_int = __builtint__.sorted(ind_lon(np.arange(1,2+1)),'ascend')
    lon_sort,iLon = __builtint__.sorted(lon_all(ind_lon_int),'ascend')
    # interpolation coefficients
    ylat = np.array([[lat_sort(1) - ipp(1)],[ipp(1) - lat_sort(2)]])
    xlon = np.array([lon_sort(2) - lon_ipp_1,lon_ipp_1 - lon_sort(1)])
    r = 1 / ((lon_sort(2) - lon_sort(1)) * (lat_sort(1) - lat_sort(2)))
    # vTEC values
    Q1 = np.array([latBlock2_t1(ind_lon_int(iLon(1))),latBlock1_t1(ind_lon_int(iLon(1))),latBlock2_t1(ind_lon_int(iLon(2))),latBlock1_t1(ind_lon_int(iLon(2)))])
    # interpolated vTEC for epoch one
    vTec_t1 = r * xlon * Q1 * ylat
    clear('ylat','xlon','ind_lon','ind_lon_int','lon_sort','iLon')
    # epoch two
# find the 2 closest longitude indices (based on the 2nd modified IPP longitude)
    __,ind_lon = __builtint__.sorted(np.abs(lon_ipp_2 - lon_all),'ascend')
    ind_lon_int = __builtint__.sorted(ind_lon(np.arange(1,2+1)),'ascend')
    lon_sort,iLon = __builtint__.sorted(lon_all(ind_lon_int),'ascend')
    # interpolation coefficients
    ylat = np.array([[lat_sort(1) - ipp(1)],[ipp(1) - lat_sort(2)]])
    xlon = np.array([lon_sort(2) - lon_ipp_2,lon_ipp_2 - lon_sort(1)])
    r = 1 / ((lon_sort(2) - lon_sort(1)) * (lat_sort(1) - lat_sort(2)))
    # vTEC values
    Q2 = np.array([latBlock2_t2(ind_lon_int(iLon(1))),latBlock1_t2(ind_lon_int(iLon(1))),latBlock2_t2(ind_lon_int(iLon(2))),latBlock1_t2(ind_lon_int(iLon(2)))])
    # interpolated vTEC for epoch two
    vTec_t2 = r * xlon * Q2 * ylat
    # interpolation in time
    vTec = ((time_sort(2) - time_dec) * vTec_t1) + ((time_dec - time_sort(1)) * vTec_t2)
    return vTec
    
    return vTec# Reza Bordbari
# This script selects the appropriate latitude block in the Rinex file

'''

def get_lat_block(IONEXFile = None,fullDate = None,time_dec = None,time_sort = None,lat_sort = None): 
    fullDate1 = fullDate
    if time_dec > 23:
        time_sort[2] = 0
        fullDate2 = fullDate + np.array([[0],[0],[1]])
    else:
        fullDate2 = fullDate
    
    # number of lines in the IONEX file.
    NoL = CalcNumOfLines(IONEXFile)
    # find the epochs and the two latitude blocks per epoch
    DataFile = open(IONEXFile,'r')
    fseek(DataFile,0,'bof')
    tline = fgetl(DataFile)
    for i in np.arange(1,NoL+1).reshape(-1):
        # epoch one
        e1 = strfind(tline,np.array([num2str(transpose(np.array([[fullDate1],[time_sort(1)],[0],[0]]))),'                        ','EPOCH OF CURRENT MAP']))
        if not (len(e1)==0) :
            condition = 1
            while condition == 1:
                # lat one
                l1 = strfind(tline,np.array([num2str(lat_sort(1),'%.1f'),'-180.0 180.0   5.0 450.0']))
                if not (len(l1)==0) :
                    a = fgetl(DataFile)
                    b = a
                    for j in np.arange(1,4+1).reshape(-1):
                        b = append(b,fgetl(DataFile))
                    latBlock1_t1 = str2num(b)
                tline = fgetl(DataFile)
                # lat two
                l2 = strfind(tline,np.array([num2str(lat_sort(2),'%.1f'),'-180.0 180.0   5.0 450.0']))
                if not (len(l2)==0) :
                    a = fgetl(DataFile)
                    b = a
                    for j in np.arange(1,4+1).reshape(-1):
                        b = append(b,fgetl(DataFile))
                    latBlock2_t1 = str2num(b)
                    condition = 2
        # epoch two
        e2 = strfind(tline,np.array([num2str(transpose(np.array([[fullDate2],[time_sort(2)],[0],[0]]))),'                        ','EPOCH OF CURRENT MAP']))
        if not (len(e2)==0) :
            condition = 1
            while condition == 1:
                # lat one
                l1 = strfind(tline,np.array([num2str(lat_sort(1),'%.1f'),'-180.0 180.0   5.0 450.0']))
                if not (len(l1)==0) :
                    a = fgetl(DataFile)
                    b = a
                    for j in np.arange(1,4+1).reshape(-1):
                        b = append(b,fgetl(DataFile))
                    latBlock1_t2 = str2num(b)
                clear('a','b','j')
                tline = fgetl(DataFile)
                # lat two
                l2 = strfind(tline,np.array([num2str(lat_sort(2),'%.1f'),'-180.0 180.0   5.0 450.0']))
                if not (len(l2)==0) :
                    a = fgetl(DataFile)
                    b = a
                    for j in np.arange(1,4+1).reshape(-1):
                        b = append(b,fgetl(DataFile))
                    latBlock2_t2 = str2num(b)
                    condition = 2
            break
        tline = fgetl(DataFile)
    return latBlock1_t1,latBlock2_t1,latBlock1_t2,latBlock2_t2


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


'''
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
'''

def calculate_daz_iono(frame, esds, framespd, method = 'gradient', out_hionos = False, out_alphas = False,
                 out_tec_master = False, out_tec_all = False, ionosource='code', use_iri_hei=False, alpha = 0.85):
    ''' Function to calculate iono correction for a given frame.

    Args:
        frame (str)                 frame ID
        esds (pandas.Dataframe)     standard esds table
        framespd (pandas.Dataframe) standard framespd table
        method (str)                gradient or liang (gradient is correct, Liang is kept from first attempts, may not work)
        out_hionos (bool)           whether to output also IRI-estimated peak height of the F2 iono layer (sadly, CODE calculates for 450 km hei which is probably wrong)
        out_tec_master (bool)
        out_tec_all (bool)
        ionosource (str)            iri or code (for IRI or GIM-based ionosphere. the latter improves RMSE!)
        use_iri_hei (bool)          if True, it estimates F2 peak altitude using IRI and uses for the correction (with any ionosource). NOTE, CODE data are for 450 km ALT..
        alpha (float or 'auto')     only for code (/jpl) source. If 'auto', it would estimate it using IRI model

    Notes: 'liang' method should include also some extra F2 height correction..
    2023/08: Liang method was first implemented here and a lot happened since that time.. Please consider it obsolete.
    '''
    # some constants
    earth_radius = 6378160 # m
    #
    if method == 'gomba': # renamed it to keep as it is
        method = 'gradient'
    if ionosource == 'iri' and (not use_iri_hei):
        print('using IRI, setting the iri to estimate F2 peak altitude')
        use_iri_hei=True
    selected_frame_esds = esds[esds['frame'] == frame].copy()
    frameta = framespd[framespd['frame']==frame]
    # extract some variables
    heading = frameta['heading'].values[0]
    scene_center_lon = frameta['center_lon'].values[0]
    scene_center_lat = frameta['center_lat'].values[0]
    # resolution_of_pixel = frameta['azimuth_resolution'].values[0]
    if 'centre_range_ok_m' in frameta:
        range_avg = frameta['centre_range_ok_m'].values[0] # 2024: GAMMA had a bug wrongly informing on centre_range (may differ by 20 km or so!). fixed most of it
    else:
        range_avg = frameta['centre_range_m'].values[0]
    master = frameta['master'].values[0]
    inc_angle_avg = frameta['avg_incidence_angle'].values[0]
    center_time=frameta['centre_time'].values[0]
    dfDC = frameta['dfDC'].values[0]
    if dfDC == 0:
        print('warning, frame '+frame+' has no info on dfDC, using default')
        dfDC = 4365 #Hz, mean value in the whole dataset
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
    # bovl_acq_dist = 7100*(2.75*2+110*0.002056) #approx. satellite velocity on the ground 7100 [m/s] * time to acquire burst overlap
    bovl_dtime = 2.75/3*2 + 0.07
    bovl_acq_dist = 7100 * bovl_dtime  # approx. satellite velocity on the ground 7100 [m/s] * time to acquire burst overlap
    # ok... we should not use PRF to calculate the time --- rather we use ratio of 110 lines of burst overlap vs 1450 lines burstwise...
    ###### do the satg_lat, lon
    azimuthDeg = heading-90 #yes, azimuth is w.r.t. N (positive to E)
    elevationDeg = 90-inc_angle_avg
    slantRange = range_avg
    #try:
    #    scene_alt = get_altitude(scene_center_lat, scene_center_lon)
    #except:
    #    scene_alt = 0
    #to get position of the satellite - UNCLEAR about slantRange - is this w.r.t. DEM? (should ask GAMMA) - if not (only in WGS-84), i should use scene_alt=0!
    # 2023/08: checked using orbits - the slantranfe is wrt ellipsoid! setting scene alt 0
    x, y, z = aer2ecef(azimuthDeg, elevationDeg, slantRange, scene_center_lat, scene_center_lon, 0) #scene_alt)
    satg_lat, satg_lon, sat_alt = ecef2latlonhei(x, y, z)
    print('debug: sat coordinates are expected as')
    print('lat,lon,alt=')
    print([satg_lat, satg_lon, sat_alt])
    print('x,y,z=')
    print([x,y,z])
    print('--------')
    Psatg = wgs84.GeoPoint(latitude=satg_lat, longitude=satg_lon, degrees=True)
    # get middle point between scene and sat - and get F2 altitude (max TEC) for it
    path = nv.GeoPath(Pscene_center.to_nvector(), Psatg.to_nvector())
    # get point in the middle
    Pmid_scene_sat = path.interpolate(0.5).to_geo_point()
    # work in dedicated table
    df = pd.DataFrame(acq_times)
    # 2024: not clear if IRI F2 peak altitude is correct. Allowing the standard 450 km assumption by CODE (I think TECs will get scaled, so the gradient should be still ok. Not tested)
    # 2025: tested. IRI is more correct. especially IRI2020
    if use_iri_hei:
        # get hionos in that middle point:
        if alpha == 'auto':
            stralpha='and alpha '
        else:
            stralpha=''
        print('extracting hmF2 '+stralpha+'estimates from IRI model')
        _, hionos, alphas = get_tecs(Pmid_scene_sat.latitude_deg, Pmid_scene_sat.longitude_deg, 800, acq_times, source='iri', returnhei = True, returnalpha=True, alpha=alpha)
        hiono_master = hionos[-1]
        selected_frame_esds['hiono'] = hionos[:-1]  ###*1000 # convert to metres, avoid last measure, as this is 'master'
        df['hiono'] = hionos
        if alpha != 'auto':
            alphas = np.array(alphas)*0 + alpha
        df['alpha'] = alphas
        # alpha_master = alphas[-1]
        selected_frame_esds['alpha'] = alphas[:-1]
    else:
        #df['hiono'] = 450 # standard altitude used by CODE
        hiono = 290  # more correct estimate
        df['hiono'] = hiono
        hionos = list(df.hiono.values)
        #hiono_master = 450
        hiono_master = hiono
        selected_frame_esds['hiono'] = hiono
        # and alphas
        if alpha == 'auto':
            print('You have set to estimate alpha but avoiding IRI - please enable use_iri_hei')
            return False
        selected_frame_esds['alpha'] = alpha
        # alpha_master = alpha
        df['alpha'] = alpha
        alphas = df.alpha.values
        #
    ############## now calculate TEC using the SLM knowledge, i.e. different A,B per epoch (!)
    # (note that the last hiono is for the master/reference epoch
    tecs_A = []
    tecs_B = []
    # do per epoch, using frame metadata (valid for reference epoch dt...)
    if 'swath_dfDC' in frameta:
        print('estimating iono gradients per swath')
        slantRange = np.array(frameta['swath_centre_range_m'].values[0])
        heading = np.array(frameta['swath_heading'].values[0])  # note: this is satellite heading, not scene heading (about 3 deg diff)
        azimuthDeg = heading - 90 #yes, azimuth is w.r.t. N (positive to E)
        theta = np.radians(np.array(frameta['swath_avg_incidence_angle'].values[0]))
        elevationDeg = 90-np.array(frameta['swath_avg_incidence_angle'].values[0])
        scene_center_lat = np.array(frameta['swath_center_lat'].values[0])
        scene_center_lon = np.array(frameta['swath_center_lon'].values[0])
        dfDC = np.array(frameta['swath_dfDC'].values[0])
        # ka = np.array(frameta['swath_ka'].values[0]) # not used....
        perswath=True
        method='gradient'
    else:
        perswath=False
    for i,a in df.iterrows():
        hiono = float(a['hiono']*1000) # m
        # overwriting input alpha param, as now fixed by pre-init above
        alpha = float(a['alpha'])
        epochdate = a['epochdate']
        #alpha_to_use = alpha
        print('running for epochtime: '+str(epochdate))
        print('assuming peak iono alt of : '+str(int(hiono/1000))+' km')
        if perswath:
            if ionosource == 'code':
                print('getting GIM VTEC')
                tecxr=get_vtec_from_code(epochdate, lat=0, lon=0, return_fullxr = True)
            range_IPP = slantRange * hiono / sat_alt
            sin_thetaiono = earth_radius/(earth_radius+hiono) * np.sin(theta)
            tecs_A_swaths = np.array([], dtype=np.float64)
            tecs_B_swaths = np.array([], dtype=np.float64)
            for j in range(len(range_IPP)):
                #print('debug: swath '+str(j))
                x, y, z = aer2ecef(azimuthDeg[j], elevationDeg[j], range_IPP[j], scene_center_lat[j], scene_center_lon[j], 0) #scene_alt)
                ippg_lat, ippg_lon, ipp_alt = ecef2latlonhei(x, y, z)
                Pippg = wgs84.GeoPoint(latitude=ippg_lat, longitude=ippg_lon, degrees=True)
                # then get A', B'
                # note, bovl_acq_dist would slightly vary per swath...
                PsatgA, _azimuth = Psatg.displace(distance=bovl_acq_dist/2, azimuth=heading[j]-180, method='ellipsoid', degrees=True)
                PsatgB, _azimuth = Psatg.displace(distance=bovl_acq_dist/2, azimuth= heading[j], method='ellipsoid', degrees=True)
                # then do intersection ... (distance is set just larger here, no meaning for bovl_acq_dist)
                PippAt, _azimuth = Pippg.displace(distance=bovl_acq_dist, azimuth=heading[j]-180, method='ellipsoid', degrees=True)
                PippBt, _azimuth = Pippg.displace(distance=bovl_acq_dist, azimuth= heading[j], method='ellipsoid', degrees=True)
                path_ipp = nv.GeoPath(PippAt, PippBt)
                Pburst_center = wgs84.GeoPoint(latitude=scene_center_lat[j], longitude=scene_center_lon[j], degrees=True)
                path_scene_satgA = nv.GeoPath(Pburst_center, PsatgA)
                path_scene_satgB = nv.GeoPath(Pburst_center, PsatgB)
                # these two points are the ones where we should get TEC
                PippA = path_ipp.intersect(path_scene_satgA).to_geo_point()
                PippB = path_ipp.intersect(path_scene_satgB).to_geo_point()
                pdist, pa1, pa2 = PippA.distance_and_azimuth(PippB, degrees=True)
                print('debug: swath '+str(j+1)+': distance between the IPP points is ' + str(int(pdist)) + ' m and their azimuth ' + str(
                    int(pa1)) + ' deg')
                if ionosource != 'code':
                    TECV_A = get_tecs(PippA.latitude_deg, PippA.longitude_deg, round(sat_alt/1000), [epochdate-pd.Timedelta(bovl_dtime/2, 's')], False, source=ionosource, alpha = alpha)[0]
                    TECV_B = get_tecs(PippB.latitude_deg, PippB.longitude_deg, round(sat_alt/1000), [epochdate+pd.Timedelta(bovl_dtime/2, 's')], False, source=ionosource, alpha = alpha)[0]
                else:
                    TECV_A = get_vtec_from_tecxr(tecxr, epochdate-pd.Timedelta(bovl_dtime/2, 's'), PippA.latitude_deg, PippA.longitude_deg, method='linear')
                    TECV_B = get_vtec_from_tecxr(tecxr, epochdate+pd.Timedelta(bovl_dtime/2, 's'), PippB.latitude_deg, PippB.longitude_deg, method='linear')
                    TECV_A = TECV_A * alpha
                    TECV_B = TECV_B * alpha
                TECS_A = TECV_A/np.sqrt(1-sin_thetaiono[j]**2)
                TECS_B = TECV_B/np.sqrt(1-sin_thetaiono[j]**2)
                tecs_A_swaths = np.append(tecs_A_swaths, TECS_A)
                tecs_B_swaths = np.append(tecs_B_swaths, TECS_B)
            tecs_A.append(tecs_A_swaths) #.mean())
            tecs_B.append(tecs_B_swaths) #.mean())
        else:
            # first, get IPP - ionosphere pierce point
            # range to IPP can be calculated using:
            # range_IPP = hiono/np.sin(theta)
            # but maybe not... let's do it simpler:
            range_IPP = slantRange * hiono / sat_alt
            x, y, z = aer2ecef(azimuthDeg, elevationDeg, range_IPP, scene_center_lat, scene_center_lon, 0) #scene_alt)
            ippg_lat, ippg_lon, ipp_alt = ecef2latlonhei(x, y, z)
            Pippg = wgs84.GeoPoint(latitude=ippg_lat, longitude=ippg_lon, degrees=True)
            # then get A', B'
            PsatgA, _azimuth = Psatg.displace(distance=bovl_acq_dist/2, azimuth=heading-180, method='ellipsoid', degrees=True)
            PsatgB, _azimuth = Psatg.displace(distance=bovl_acq_dist/2, azimuth= heading, method='ellipsoid', degrees=True)
            # then do intersection ... (distance is set just larger here, no meaning for bovl_acq_dist)
            PippAt, _azimuth = Pippg.displace(distance=bovl_acq_dist, azimuth=heading-180, method='ellipsoid', degrees=True)
            PippBt, _azimuth = Pippg.displace(distance=bovl_acq_dist, azimuth= heading, method='ellipsoid', degrees=True)
            path_ipp = nv.GeoPath(PippAt, PippBt)
            path_scene_satgA = nv.GeoPath(Pscene_center, PsatgA)
            path_scene_satgB = nv.GeoPath(Pscene_center, PsatgB)
            # these two points are the ones where we should get TEC
            PippA = path_ipp.intersect(path_scene_satgA).to_geo_point()
            PippB = path_ipp.intersect(path_scene_satgB).to_geo_point()
            pdist, pa1, pa2 = PippA.distance_and_azimuth(PippB, degrees=True)
            print('debug: distance between the IPP points is '+str(int(pdist))+' m and their azimuth '+str(int(pa1))+' deg')
            ######### get TECS for A, B
            TECV_A = get_tecs(PippA.latitude_deg, PippA.longitude_deg, round(sat_alt/1000), [epochdate], False, source=ionosource, alpha = alpha)[0]
            TECV_B = get_tecs(PippB.latitude_deg, PippB.longitude_deg, round(sat_alt/1000), [epochdate], False, source=ionosource, alpha = alpha)[0]
            # get inc angle at IPP - see iono. single layer model function
            sin_thetaiono = earth_radius/(earth_radius+hiono) * np.sin(theta)
            TECS_A = TECV_A/np.sqrt(1-sin_thetaiono**2)
            TECS_B = TECV_B/np.sqrt(1-sin_thetaiono**2)
            tecs_A.append(TECS_A)
            tecs_B.append(TECS_B)
    #
    tec_A_master = tecs_A[-1]
    tec_B_master = tecs_B[-1]
    print('debug: tec_A_master are:')
    print(tec_A_master)
    tecs_A = tecs_A[:-1]
    tecs_B = tecs_B[:-1]
    #
    #selected_frame_esds['TECS_A'] = tecs_A
    #selected_frame_esds['TECS_B'] = tecs_B
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
    #
    if perswath:
        nobursts = frame.split('_')[2]
        nobursts = [int(nobursts[:2]), int(nobursts[2:4]), int(nobursts[4:6])]
        while 0 in nobursts:
            nobursts.remove(0)
        nobursts = np.array(nobursts)
        nobovls = nobursts-1
        tecs_A2tb = []
        tecs_B2tb = []
        for t in tecs_A:
            tecs_A2tb.append(t.mean())
        for t in tecs_B:
            tecs_B2tb.append(t.mean())
        selected_frame_esds['TECS_A'] = tecs_A2tb
        selected_frame_esds['TECS_B'] = tecs_B2tb
        selected_frame_esds['daz_iono'] = 0.0
        daz_iono = []
        #
        for tecind in range(len(tecs_A)):
            tecs_A_ep = tecs_A[tecind]
            tecs_B_ep = tecs_B[tecind]
            #print('debug, tecs_A_ep:')
            #print(tecs_A_ep)
            daz_iono_ep = []
            for j in range(len(tec_A_master)):
                tecovl = (tec_A_master[j] - tecs_A_ep[j])/fH[j] - (tec_B_master[j] - tecs_B_ep[j])/fL[j]
                daz_iono_sw = 2*PRF*k/c/dfDC[j] * tecovl
                daz_iono_ep.append(daz_iono_sw*nobovls[j])
            # now get it as avg for no of bovls:
            daz_iono.append(np.array(daz_iono_ep).sum()/nobovls.sum())
        selected_frame_esds['daz_iono'] = np.array(daz_iono)
        daz_iono = selected_frame_esds['daz_iono']
        # return mean values for later use
        #tec_A_master = np.mean(tec_A_master)
        #tec_B_master = np.mean(tec_B_master)
        # np.array(tecs_A).mean(axis=1)
        tec_A_master = (tec_A_master*nobovls).sum()/nobovls.sum()
        tec_B_master = (tec_B_master*nobovls).sum()/nobovls.sum()
        tecs_A = (tecs_A*nobovls).sum(axis=1)/nobovls.sum()
        tecs_B = (tecs_B*nobovls).sum(axis=1)/nobovls.sum()
    else:
        selected_frame_esds['TECS_A'] = tecs_A
        selected_frame_esds['TECS_B'] = tecs_B
        if method == 'gradient':
            # 08/2021 - empirically checked, correct:
            #tecovl = (selected_frame_esds['TECS_B'] - tec_B_master)/(fL*fL) - (selected_frame_esds['TECS_A'] - tec_A_master)/(fH*fH)
            #daz_iono = 2*PRF*k*f0/c/dfDC * tecovl
            # 04/2022 - actually the squares seem not needed, based directly on iono2phase (see article):
            # note it is 'master' minus 'epoch S', since i start from ifg as M * complex conjugate of S (so basically phase M - phase S, phase is negative for higher delay)
            # 2023 - oh well, it was best to do straightforward equations, see Lazecky et al., 2023, GRL
            tecovl = (tec_A_master - selected_frame_esds['TECS_A'])/fH - (tec_B_master - selected_frame_esds['TECS_B'])/fL
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
            # or relate it to burst overlap interval?
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
    out = [daz_iono]
    if out_hionos:
        out.append(hionos)
    if out_alphas:
        out.append(alphas)
    if out_tec_all or out_tec_master:
        out.append(tec_A_master)
        out.append(tec_B_master)
    if out_tec_all:
        out.append(tecs_A)
        out.append(tecs_B)
    if len(out)==1:
        out=out[0]
    return out



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
