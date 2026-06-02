#!/usr/bin/env python3
'''
import os, glob
import daz_lib_licsar as dll
for tr in range(175):
    trfrs=glob.glob(str(tr+1)+'/???[A,D]_?????_??????')
    for trfr in trfrs:
        fr=trfr.split('/')[1]
        print(fr)
        outcsv = os.path.join(trfr, 'metadata', fr+'.azirg.csv')
        if not os.path.exists(outcsv):
            try:
                a=dll.get_daz_frame(fr, True, True)
                a[['epoch','daz','cc_azi','cc_range','daz_iono','daz_SET','drg_iono_mm', 'dTECS', 'drg_SET_mm']].to_csv(outcsv)
            except:
                print('error with frame '+fr)

'''
import LiCSquery as lq
import numpy as np
from LiCSAR_lib.LiCSAR_misc import *
import os, glob
import pandas as pd
import framecare as fc
try:
    import rioxarray
except:
    print('no rioxarray module - GACOS TIF files will not be read (for rng correction)')

try:
    from orbit_lib import *
except:
    print('LiCSAR orbit_lib not found, cannot process orbit files')


from daz_lib import *
import datetime as dt # just in case..

def extract_all2txt(outfr = 'frames.txt', outdaz = 'esds.txt', inframelist = None, fix_epoch_time = False):
    """ Main function to extract all frame and daz data from the LiCSAR database.

    If inframelist is None, it will import all existing (initialised) frames. Otherwise you can just add a list of frames, e.g.
    inframelist = ['117A_04827_031113']
    """
    frames = create_framelist(outfile = outfr, inframelist = inframelist)
    extract2txt_esds_all_frames(framelist = frames, outfile=outdaz, fix_epoch_time = fix_epoch_time)
    print('done, please continue by daz_01_prepare_inputs.py')


def create_framelist(outfile='frames.txt', inframelist = None):
    """ Creates frames.txt for all frames in LiCSInfo that contains table in:
    frame,master,center_lon,center_lat

    If inframelist is None, it will import all existing (initialised) frames. Otherwise you can just add a list of frames, e.g.
    inframelist = ['117A_04827_031113']
    """
    if type(inframelist) == type(None):
        print('ingesting all frames')
        framespd = fc.get_all_frames(only_initialised = True, merge = True)
    else:
        print('acquiring selected frames info')
        for fr in inframelist:
            framespd = fc.get_frames_gpd(inframelist)
    # get lons,lats:
    c=framespd.geometry.centroid
    lons = []
    lats = []
    for cc in c:
        lons.append(cc.coords[0][0])
        lats.append(cc.coords[0][1])
    # get masters:
    masters = []
    for frame in framespd.frameID:
        masters.append(fc.get_master(frame))
    framespd['frame'] = framespd['frameID']
    framespd['master'] = masters
    framespd['center_lon'] = lons
    framespd['center_lat'] = lats
    framespd = framespd.drop(columns=['geometry', 'frameID'])
    framespd.to_csv(outfile, index=False)
    return list(framespd.frame)

'''
e.g.
framelist=pd.read_csv('frames.txt'); framelist=list(framelist.frame)

where you may get frames.txt using the function create_framelist(), or using the previous code:
cd $LiCSAR_procdir
for tr in `seq 1 175`; do for f in `ls $tr`; do
  m=`ls $tr/$f/SLC | head -n1`; 
  hgtfile=$LiCSAR_public/$tr/$f/metadata/$f'.geo.hgt.tif'; 
  ll=`gdalinfo $hgtfile | grep ^Center`; 
  lon=`echo $ll | cut -d "," -f1 | cut -d '(' -f2`; 
  lat=`echo $ll | cut -d ")" -f1 | cut -d ',' -f2`; 
  echo $f","$m","$lon","$lat >> $outfr;  
done;done
'''

def extract2txt_esds_all_frames(framelist, outfile='esds.txt', fix_epoch_time = False):
    dazes=pd.DataFrame()
    for frame in framelist:
        try:
            a=extract2txt_esds_frame(frame, fix_epoch_time = fix_epoch_time)
            #dazes=dazes.append(a)
            dazes = pd.concat([dazes, a], ignore_index=True)
        except:
            print('frame '+frame+' is empty')
    dazes=dazes.reset_index(drop=True)
    dazes.to_csv(outfile, index=False)


def extract2txt_esds_frame(frame, fix_epoch_time = False, datemin=dt.date(2014,10,1), datemax=dt.date.today()):
    '''
    extracts to esds txt full data for given frame, from database
    the resultant txt file is a csv as:
    frame,esd_master,epoch,daz_total_wrt_orbits,daz_cc_wrt_orbits,orbits_precision,version

    or with epochtime as well, if this was set by fix_epoch_time (takes long...)
    '''
    a = get_daz_frame(frame, datemin=datemin, datemax=datemax)
    if 'epoch' in a:
        a['epoch']=a.epoch.apply(lambda x: x.strftime('%Y%m%d'))
    a['esd_master']=a.rslc3.apply(lambda x: x.strftime('%Y%m%d'))
    a['daz_total_wrt_orbits']=a.daz  #+a.cc_azi daz is the final value (ICC+SD)
    a['orbits_precision'] = 'P'  # only Ps should be in database
    a['version'] = 'm' # i forgot what this is for, but should be ok any letter (?)
    a['frame'] = frame
    a=a.rename(columns={'cc_azi':'daz_cc_wrt_orbits'})
    cols = ['frame','esd_master','epoch','daz_total_wrt_orbits','daz_cc_wrt_orbits','orbits_precision','version']
    if fix_epoch_time:
        bperp, a['epochtime'] = fc.estimate_bperps(frame, list(a.epoch.values), return_epochsdt = True)
        cols = cols + ['epochtime']
    return a[cols]


def get_platemotion_en(df, collon = 'centroid_lon', collat = 'centroid_lat', outcolnm='eur', plate = 'Eurasia'):
    import licsbas_mintpy_PMM as pmm
    #
    # getting plate data
    plate = pmm.ITRF2014_PMM[plate]
    pole_obj = pmm.EulerPole(
        wx=plate.omega_x,
        wy=plate.omega_y,
        wz=plate.omega_z,
        unit='mas/yr',
    )
    #
    lats = np.array(df[collat].values)
    lons = np.array(df[collon].values)
    # finally getting the plate velocities from the Euler pole definition over the frame area
    ve, vn, vu = pole_obj.get_velocity_enu(lats, lons, alt=0.0, ellps=True)
    df[outcolnm+'_E'] = ve*1000
    df[outcolnm + '_N'] = vn*1000 # in mm/y
    return df


def estimate_vels_frame(frame, tocsvs = True):
    ''' this will perform full processing... if frame.csv is present in the folder, it will load that one instead of full regen'''
    csvfile = frame+'.csv'
    if os.path.exists(csvfile):
        esds = pd.read_csv(csvfile)
        try:
            esds = esds.drop('Unnamed: 0', axis=1)
        except:
            pass
    else:
        esds=get_daz_frame(frame, fulloutput=True, include_corrections=True)
    esds = esds[esds['cc_range']!=0] # really needed?
    esds = esds[esds['daz'] != 0]  # really needed?
    esds['epochdate'] = esds['epoch'].copy(deep=True)
    esds = esds.drop('epoch', axis=1)
    esds['epochdate'] = esds.apply(lambda x : pd.to_datetime(str(x.epochdate)).date(), axis=1)
    esds = esds[esds['epochdate'] > dt.date(2016,3,1)]
    esds['years_since_beginning'] = 0.0
    firstdatei = esds['epochdate'].min()
    esds['years_since_beginning'] = esds['epochdate'] - firstdatei
    esds['years_since_beginning'] = esds['years_since_beginning'].apply(lambda x: float(x.days)/365.25)
    frameta=get_frameta(frame)
    esds['drg_mm_notide_noiono_nogacos']=esds['cc_range']*float(frameta['range_resolution']*1000) - esds['drg_SET_mm'] - esds['drg_iono_mm'] - esds['drg_GACOS_mm']
    esds['daz_mm_notide_noiono']=(esds['daz']-esds['daz_iono']-esds['daz_SET'])*float(frameta['azimuth_resolution']*1000)
    mdate=pd.Timestamp(frameta.master.values[0])
    esds['S1AorB'] = flag_s1b(esds.epochdate.values, mdate, 'A', True)
    import daz_timeseries as dts
    esds['drg_final_mm'] = esds['drg_mm_notide_noiono_nogacos'].copy(deep=True)
    esds['daz_final_mm'] = esds['daz_mm_notide_noiono'].copy(deep=True)
    esdsB = esds[esds['S1AorB']=='B']
    v,c,stderr,c_AB = dts.estimate_s1ab(esds, col = 'daz_mm_notide_noiono', rmsiter = 50, printout = True)
    frameta['vel_az_mmy'] = [v]
    frameta['intercept_az'] = [c]
    frameta['vel_az_stderr_mmy'] = [stderr]
    frameta['S1AB_offset_az'] = [c_AB]
    esdsB['daz_final_mm'] = esdsB['daz_final_mm'] - c_AB
    #
    v,c,stderr,c_AB = dts.estimate_s1ab(esds, col = 'drg_mm_notide_noiono_nogacos', rmsiter = 50, printout = True)
    frameta['vel_rg_mmy'] = [v]
    frameta['intercept_rg'] = [c]
    frameta['vel_rg_stderr_mmy'] = [stderr]
    frameta['S1AB_offset_rg'] = [c_AB]
    esdsB['drg_final_mm'] = esdsB['drg_final_mm'] - c_AB
    esds.update(esdsB)
    # finally re-estimate vels
    selesds = esds[np.abs(esds['daz_final_mm']-esds['daz_final_mm'].median())<400]
    x = selesds.years_since_beginning.values  # .transpose()
    X = np.array([x]).transpose()
    y = selesds['daz_final_mm'].values
    huber = HuberRegressor(alpha = 1, epsilon = 1.35)
    huber.fit(X, y)
    frameta['vel_az_mmy_huber'] = [huber.coef_[0]]
    #
    selesds = esds[np.abs(esds['drg_final_mm'] - esds['drg_final_mm'].median()) < 400]
    x = selesds.years_since_beginning.values #.transpose()
    X = np.array([x]).transpose()
    y = selesds['drg_final_mm'].values
    huber = HuberRegressor(alpha=1, epsilon=1.35)
    huber.fit(X, y)
    frameta['vel_rg_mmy_huber'] = [huber.coef_[0]]
    if tocsvs:
        esds.to_csv(frame+'.esds.csv')
        frameta.to_csv(frame+'.vels.csv')
    return esds, frameta


def get_daz_frame(frame, fulloutput = True, include_corrections = False, use_iri_hei = False, corr_per_swath = True,
                  datemin=dt.date(2014,10,1), datemax=dt.date.today()):
    ''' Function to extract all frame daz values from the LiCSInfo database.

    Args:
        frame (str)                 LiCSAR frame ID
        fulloutput (bool)           if True, will return all information from the database (such as cc_rg), otherwise only daz values [mm]
        include_corrections (bool)  if True, will perform also SET and iono corrections and add to the table (or daz if False fulloutput)
        use_iri_hei (bool)          only with include_corrections - will try scale TEC values with IRI model
        corr_per_swath (bool)       will bit improve corrections using center region per swath rather than from frame center
        datemin, datemax (dt.date)  will limit the epochs
    '''
    polyid=lq.get_frame_polyid(frame)[0][0]
    daztb = lq.do_pd_query('select * from esd where polyid={};'.format(polyid))
    daztb = daztb[daztb.epoch >= datemin]
    daztb = daztb[daztb.epoch <= datemax]
    if include_corrections:
        frameta = get_frameta(frame, perswath=corr_per_swath)
        mastr = str(pd.to_datetime(frameta.master[0]).date())
        esds = extract2txt_esds_frame(frame, datemin=datemin, datemax=datemax)
        esds['epochdate'] = esds.apply(lambda x : pd.to_datetime(str(x.epoch)).date(), axis=1)
        import daz_iono as di
        if use_iri_hei:
            alpha = 'auto'
        else:
            alpha = 0.85
        daz_iono, hionos, alphas, tec_A_master, tec_B_master, tecs_A, tecs_B = di.calculate_daz_iono(frame, esds, frameta, method='gradient', out_hionos=True, out_alphas=True,
                                                                                     out_tec_all=True, ionosource='code', use_iri_hei=use_iri_hei, alpha = alpha)
        daztb['daz_iono'] = daz_iono
        # SET will still work frame-centered only (need this to change??)
        print('Getting solid earth tides corrections')
        daztb['daz_SET'] = get_SET_for_frame_dazes(frameta, esds)
        if fulloutput:
            import rioxarray
            # getting drg iono:
            tecs = (np.array(tecs_A) + np.array(tecs_B)) / 2
            k = 40.308193  # m^3 / s^2
            f0 = 5.4050005e9  # 1/s
            tec_ref = (tec_A_master + tec_B_master) / 2
            daztb['drg_iono_mm'] = 1000* (tecs - tec_ref) * k / (f0 * f0) # in mm
            daztb['dTECS'] = tecs - tec_ref  # for later correlation to azi shift (seems abs TEC plays role!)
            daztb['alpha'] = alphas[:-1]
            daztb['hmF2'] = hionos[:-1]
            # SET for drg --- we need to get ENUs:
            geoframedir = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])), frame)
            e = os.path.join(geoframedir, 'metadata', frame + '.geo.E.tif')
            n = os.path.join(geoframedir, 'metadata', frame+'.geo.N.tif')
            u = os.path.join(geoframedir, 'metadata', frame + '.geo.U.tif')
            lon = frameta['center_lon'][0]
            lat = frameta['center_lat'][0]
            centre_time = frameta['centre_time'][0]
            try:
                e = rioxarray.open_rasterio(e).squeeze('band').drop('band')
                n = rioxarray.open_rasterio(n).squeeze('band').drop('band')
                u = rioxarray.open_rasterio(u).squeeze('band').drop('band')
                e = float(e.where(e != 0).mean())
                n = float(n.where(n != 0).mean())
                u = float(u.where(u != 0).mean())
                drg_tides_mm = []
                gcs_mm = []
                for edt in esds.epochdate.values:
                    epochdtstr = str(edt) + 'T' + centre_time
                    E, N, U = get_SET_coords(lon, lat, epochdtstr)
                    drg_tide_m = E*e + N*n + U*u
                    drg_tides_mm.append(1000*drg_tide_m*(-1))
                    gcs = get_gacos_in_coord(lon,lat,str(edt).replace('-',''), frame, inmm=True)
                    # gcs = get_gacos_in_coord(lon, lat, str(edt).replace('-',''), frame, inmm=True, domean=True) # mean is WORSE!!!!
                    gcs_mm.append(gcs*(-1))
                drg_tides_mm = np.array(drg_tides_mm)
                gcs_mm = np.array(gcs_mm)
                daztb['drg_SET_mm'] = drg_tides_mm
                daztb['drg_GACOS_mm'] = gcs_mm
                try:
                    mastrT = mastr + 'T' + centre_time
                    E, N, U = get_SET_coords(lon, lat, mastrT)
                    mastr_SET = 1000*(E * e + N * n + U * u * (-1))
                    mastr_GACOS = get_gacos_in_coord(lon, lat, mastr.replace('-',''), frame, inmm=True, domean=False)
                    mastr_GACOS = mastr_GACOS * (-1)
                    daztb['drg_SET_mm'] = daztb['drg_SET_mm'] - mastr_SET
                    daztb['drg_GACOS_mm'] = daztb['drg_GACOS_mm'] - mastr_GACOS
                except:
                    print('correction wrt reference epoch did not succeed - expect overall offset')
            except:
                print('some issue getting range SET corrections..')
            return daztb
        else:
            # if we want to return only corrected daz_mm
            daz_mm = daztb.set_index('epoch')['daz'] * 14000
            daz_mm.values = daz_mm.values - daztb['daz_iono'].values * 14000 - daztb['daz_SET'].values * 14000
            return daz_mm
    else:
        if fulloutput:
            return daztb
        else:
            daz_mm = daztb.set_index('epoch')['daz'] * 14000
            return daz_mm



def get_center_vel(parfile):
    center_time=get_param_gamma('center_time', parfile, floatt = True, pos = 0)
    if not os.path.exists(parfile+'.orb'):
        rc = os.system("ORB_prop_SLC "+parfile+" - - - 1 >/dev/null; ORB_prop_SLC "+parfile+" - - - 1 | grep 'output sv' > "+parfile+".orb")
    #time.sleep(0.5)
    sv = pd.read_csv(parfile+'.orb', delim_whitespace=True,header=None)
    svvs = []
    for i in sv[6]:
        svvi=get_param_gamma('state_vector_velocity_'+str(i), parfile, pos=0)
        #svvi=get_param_gamma('state_vector_velocity_'+str(i), parfile, pos=2)
        svvs.append(svvi)
    sv['vel1']=svvs
    # least squares to get vel for given center_time:
    x=sv[3].values
    A = np.vstack([x, np.ones(len(x))]).T
    y = sv['vel1'].values
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]
    center_vel=m*center_time+c
    return center_vel


def get_velocities_per_sat(rslcdir='RSLC'):
    epochs=os.listdir(rslcdir)
    epochspd=pd.DataFrame(epochs)
    sats=[]
    vel1s=[]
    for epoch in epochs:
        print(epoch)
        parfile=os.path.join(rslcdir,epoch,epoch+'.rslc.par')
        if os.path.exists(parfile):
            sat=get_param_gamma('sensor', parfile, floatt = False, pos = 0)
            vel1=get_center_vel(parfile)
        else:
            sat=''
            vel1=np.nan
        sats.append(sat)
        vel1s.append(vel1)
    epochspd['sat']=sats
    epochspd['vel1']=vel1s
    epochspd=epochspd.dropna()
    epochspd['epochdate'] = epochspd.apply(lambda x : pd.to_datetime(str(x[0])).date(), axis=1)
    epochspd=epochspd.set_index(epochspd['epochdate']).sort_index()
    return epochspd


'''
import matplotlib.pyplot as plt
vel2=epochspd.copy(deep=True)
epochspd=get_velocities_per_sat(rslcdir)
pp=epochspd[epochspd['sat']=='S1A'].plot()
epochspd[epochspd['sat']=='S1B'].plot(ax=pp)
plt.show()
'''



'''
how to get s1a/b? try this: BUT CAREFUL! daz might be just daztb.daz (not plus cc_azi - check it first)
frame = 
polyid=fc.lq.get_frame_polyid(frame)[0][0]
daztb = fc.lq.do_pd_query('select * from esd where polyid={};'.format(polyid))
daztb = daztb.set_index(daztb.epoch).sort_index()
daz = (daztb['daz'] + daztb['cc_azi'])*14000

daz=daz[daz<1000][daz>-1000]
s1bs = []
s1as = []
for e in daz.index:
    if (np.mod((e - daz.index[0]).days, 12) == 0):
        s1as.append(e)
    else:
        s1bs.append(e)
B = daz[np.isin(daz.index, s1bs)]
A = daz[np.isin(daz.index, s1as)]

'''
def get_azshift_SD(offile):
    azshift_SD = float(grep1line('azimuth_offset_polynomial', offile).split()[1])
    return azshift_SD

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

# orig function - now in framecare
def get_frame_master_s1ab(frame):
    tr = int(frame[:3])
    metafile = os.path.join(os.environ['LiCSAR_public'], str(tr), frame, 'metadata', 'metadata.txt')
    if not os.path.exists(metafile):
        print('metadata file does not exist for frame '+frame)
        return 'X'
    try:
        primepoch = grep1line('master=',metafile).split('=')[1]
    except:
        print('metadata file does not contain prim epoch info for frame '+frame)
        return 'X'
    path_to_slcdir = os.path.join(os.environ['LiCSAR_procdir'], str(tr), frame, 'SLC', primepoch)
    try:
        out = os.path.basename(glob.glob(path_to_slcdir+'/S1*')[0])[2]
    except:
        print('error getting the value for frame '+frame)
        out = 'X'
    return out


def extract_frame_master_s1abs(framespd):
    s1abs = []
    for frame in framespd.frame:
        s1abs.append(get_frame_master_s1ab(frame))
    framespd['S1AorB'] = s1abs
    return framespd


def get_frameta(frame, perswath=False):
    a = pd.DataFrame()
    tr = int(frame[:3])
    metafile = os.path.join(os.environ['LiCSAR_public'], str(tr), frame, 'metadata', 'metadata.txt')
    if not os.path.exists(metafile):
        print('metadata file does not exist for frame ' + frame)
        return False
    primepoch = grep1line('master=', metafile).split('=')[1]
    path_to_slcdir = os.path.join(os.environ['LiCSAR_procdir'], str(tr), frame, 'SLC', primepoch)
    try:
        dfDC, ka = get_dfDC(path_to_slcdir, returnperswath = perswath)
        if perswath:
            a['swath_ka'] = [ka]
            a['swath_dfDC'] = [dfDC]
            dfDC, ka = get_dfDC(path_to_slcdir) # this is to keep per-frame values..
    except:
        print('some error occurred during extracting dfDC per swath - cancelling per swath option')
        perswath = False
        try:
            dfDC, ka = get_dfDC(path_to_slcdir)
        except:
            print('some error occurred getting average dfDC for frame ' + frame+' - setting default dfDC and ka')
            dfDC = 4365
            ka = -2000
    a['frame'] = [frame]
    a['master'] = [primepoch]
    a['ka'] = [ka]
    a['dfDC'] = [dfDC]
    if perswath:
        parfiles = glob.glob(os.path.join(path_to_slcdir,primepoch)+'.IW?.slc.par')
        heading = []
        azimuth_resolution = []
        avg_incidence_angle = []
        centre_range_m = []
        centre_time = []
        lon = []
        lat = []
        for par in parfiles:
            heading.append(get_param_gamma('heading',par))
            azimuth_resolution.append(get_param_gamma('azimuth_pixel_spacing',par))
            avg_incidence_angle.append(get_param_gamma('incidence_angle',par))
            # old GAMMA had issue with centre_range - just to be sure, averaging from near and far range
            nearrange = get_param_gamma('near_range_slc',par)
            farrange = get_param_gamma('far_range_slc',par)
            centre_range_m.append((nearrange+farrange)/2)
            cdate = grep1line('date', par).split()[4:]
            cdate = cdate[0]+':'+cdate[1]+':'+cdate[2]
            centre_time.append(cdate)
            lon.append(get_param_gamma('center_longitude',par))
            lat.append(get_param_gamma('center_latitude',par))
        a['swath_heading'] = [heading]
        a['swath_azimuth_resolution'] = [azimuth_resolution]
        a['swath_avg_incidence_angle'] = [avg_incidence_angle]
        a['swath_centre_range_m'] = [centre_range_m]
        a['swath_centre_time'] = [centre_time]
        a['swath_center_lon'] = [lon]
        a['swath_center_lat'] = [lat]
    heading = float(grep1line('heading', metafile).split('=')[1])
    try:
        azimuth_resolution = float(grep1line('azimuth_resolution', metafile).split('=')[1])
        range_resolution = float(grep1line('range_resolution', metafile).split('=')[1])
    except:
        print('WARNING, no resolution info in metafile '+metafile+' - using defaults')
        azimuth_resolution = 14.0 #m
        range_resolution = 2.3
    avg_incidence_angle = float(grep1line('avg_incidence_angle', metafile).split('=')[1])
    try:
        centre_range_m = float(grep1line('centre_range_ok_m', metafile).split('=')[1])
    except:
        centre_range_m = float(grep1line('centre_range_m', metafile).split('=')[1])
    centre_time = grep1line('center_time', metafile).split('=')[1]
    framegeom = fc.get_frames_gpd([frame])
    c = framegeom.geometry.centroid
    lon = c[0].coords[0][0]
    lat = c[0].coords[0][1]
    a['heading'] = [heading]
    a['azimuth_resolution'] = [azimuth_resolution]
    a['range_resolution'] = [range_resolution]
    a['avg_incidence_angle'] = [avg_incidence_angle]
    a['centre_range_m'] = [centre_range_m]
    a['centre_time'] = [centre_time]
    a['center_lon'] = [lon]
    a['center_lat'] = [lat]
    try:
        hei = grep1line('avg_height', metafile).split('=')[1]
    except:
        print('no height information, returning 0 for frame ' + frame)
        hei = 0
    a['avg_height'] = [hei]
    return a


def generate_framespd(fname = 'esds2021_frames.txt', outcsv = 'framespd_2021.csv'):
    ''' Function to collect additional data for frames listed in fname txt file, and store as a csv.

    Note: input fname is generated using create_framelist and has header:
    frame, master, center_lon, center_lat

    Output - csv with header:
    frame, master, center_lon, center_lat, heading, azimuth_resolution, avg_incidence_angle, centre_range_m, centre_time, ka, dfDC, avg_height, S1AorB
    '''
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
    #a['kr']=0.00
    a['dfDC'] = 0.00
    a['avg_height'] = 0.00
    for i,row in a.iterrows():
        if np.mod(i, 100)==0:
            print('Processed '+str(i)+'/'+str(len(a))+' frames', flush=True)
        frame=row['frame']
        #print(frame)
        tr = int(frame[:3])
        metafile = os.path.join(os.environ['LiCSAR_public'], str(tr), frame, 'metadata', 'metadata.txt')
        if not os.path.exists(metafile):
            print('metadata file does not exist for frame '+frame)
            continue
        try:
            primepoch = grep1line('master=',metafile).split('=')[1]
        except:
            print('the frame '+frame+' has no information on primary epoch in metadata.txt. Skipping')
            continue
        path_to_slcdir = os.path.join(os.environ['LiCSAR_procdir'], str(tr), frame, 'SLC', primepoch)
        # 
        #if frame == '174A_05407_121212':
        #    heading = -10.157417
        #    azimuth_resolution = 13.968690
        #    avg_incidence_angle = 39.5118
        #    centre_range_m = 878941.4133
        #    centre_time = '14:52:00'
    #   #     kt = 
        try:
            heading = float(grep1line('heading',metafile).split('=')[1])
            azimuth_resolution = float(grep1line('azimuth_resolution',metafile).split('=')[1])
            avg_incidence_angle = float(grep1line('avg_incidence_angle',metafile).split('=')[1])
            try:
                centre_range_m = float(grep1line('centre_range_ok_m',metafile).split('=')[1])
            except:
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
            #dfDC, ka, kr = get_dfDC(path_to_slcdir)
            dfDC, ka = get_dfDC(path_to_slcdir)
        except:
            print('some error occurred during frame '+frame)
            dfDC = 0
            ka = 0
            #kr = 0
        try:
            hei = grep1line('avg_height',metafile).split('=')[1]
        except:
            print('no height information, returning 0 for frame '+frame)
            hei = 0
        a.at[i,'heading'] = heading
        a.at[i,'azimuth_resolution']  = azimuth_resolution
        a.at[i,'avg_incidence_angle']  = avg_incidence_angle
        a.at[i,'centre_range_m']  = centre_range_m
        a.at[i,'centre_time']  = centre_time
    #    a.at[i,'kt']  = kt
        a.at[i,'dfDC']  = dfDC
        a.at[i,'ka']  = ka
        a.at[i,'avg_height'] = hei
        #a.at[i,'kr']  = kr
    print('Information on frames collected. Flagging satellite ID (S1A/B) of the reference epoch')
    a = extract_frame_master_s1abs(a)
    print('cleaning the dataset')
    preclean = len(a)
    try:
        a = clean_framespd(a)
    except:
        print('some error during cleaning, keeping the original table')
    postclean = len(a)
    print('{0}/{1} frames remain after the cleaning'.format(str(postclean), str(preclean)))
    a.to_csv(outcsv, float_format='%.4f', index=False)
    return a


def clean_framespd(a):
    a=a[a['S1AorB']!='X']
    a=a[a['master']!='False']
    a = a[a['heading'] != 0]
    a = a[a['azimuth_resolution'] != 0]
    a = a[a['avg_incidence_angle'] != 0]
    a = a[a['centre_range_m'] != 0]
    a = a[a['centre_time'] != 0]
    a = a[a['ka'] != 0]
    a = a[a['dfDC'] != 0]
    a = a[a['frame'].str[5]!='S']
    for fie in a:
        a = a[~np.isnan(a[fie])]
    a = a.reset_index(drop=True)
    return a


def get_avg_height(frame):
    tr = int(frame[:3])
    metafile = os.path.join(os.environ['LiCSAR_public'], str(tr), frame, 'metadata', 'metadata.txt')
    if not os.path.exists(metafile):
        print('metadata file does not exist for frame '+frame)
        return np.nan
    try:
        hei = grep1line('avg_height=',metafile).split('=')[1]
    except:
        print('no info on height for frame '+frame)
        return np.nan
    return hei


def get_avg_height_framespd(framespd):
    aaa = []
    for i,r in framespd.iterrows():
        frame=r['frame']
        try:
            val=float(get_avg_height(frame))
        except:
            val = np.nan
        aaa.append(val)
    framespd['avg_height'] = aaa
    return framespd


def get_dfDC(path_to_slcdir, f0=5405000500, burst_interval = 2.758277, returnka = True, returnperswath = False):
    #f0 = get_param_gamma('radar_frequency', parfile)
    #burst_interval = get_param_gamma('burst_interval', topsparfile)
    epoch = os.path.basename(path_to_slcdir)
    frame = path_to_slcdir.split('/')[-3]
    # parfile = os.path.join(path_to_slcdir, epoch+'.slc.par')
    #parfile = glob.glob(path_to_slcdir+'/????????.slc.par')[0]
    #topsparfiles = glob.glob(path_to_slcdir+'/????????.IW?.slc.TOPS_par')
    #iwparfiles = glob.glob(path_to_slcdir+'/????????.IW?.slc.par')
    #
    lam = speed_of_light / f0
    dfDC = []
    kas = []
    numbursts = [ int(frame.split('_')[2][:2]), int(frame.split('_')[2][2:4]), int(frame.split('_')[2][4:6])]
    #krs = []
    #print('This is a proper solution but applied to primary SLC image. originally it is applied by GAMMA on the RSLC...')
    #for n in range(len(topsparfiles)):
    for n in [1,2,3]:
        topsparfile = os.path.join(path_to_slcdir, epoch+'.IW'+str(n)+'.slc.TOPS_par')
        iwparfile = os.path.join(path_to_slcdir, epoch+'.IW'+str(n)+'.slc.par')
        if (not os.path.exists(iwparfile)) or (not os.path.exists(topsparfile)):
            dfDC.append(np.nan)
            kas.append(np.nan)
            numbursts[n-1] = np.nan
        else:
            #topsparfile = topsparfiles[n]
            #iwparfile = iwparfiles[n]
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
    if not returnperswath:
        numbursts = np.array(numbursts)
        dfDC = np.nansum(numbursts*np.array(dfDC)) / np.sum(numbursts)
        kas = np.nansum(numbursts*np.array(kas)) / np.sum(numbursts)
    #kr = np.mean(krs)
    if returnka:
        return dfDC, kas #, kr
    else:
        return dfDC


# to include rdc_trans results (e.g. 099A comparison figure in the article):

def get_shifts_from_qualfile(qualfile):
    daz_icc = 0
    dr_icc = 0
    gr = grep_full('matching_iteration_',qualfile)
    if gr:
        for l in gr:
            if 'daz' in l:
                daz_icc = daz_icc + float(l.split()[1])
                dr_icc = dr_icc + float(l.split()[2])
        daz_sd = float(grep1line('Total azimuth offset', qualfile).split(':')[1].split()[0])
    else:
        # different version...
        gr = grep_full('^daz =',qualfile)
        if gr:
            for l in gr:
                daz_icc = daz_icc + float(l.split()[2])
            grr = grep_full('^dr =',qualfile)
            for l in grr:
                dr_icc = dr_icc + float(l.split()[2])
            daz_sd = float(grep1line('Total azimuth offset', qualfile).split(':')[1].split()[0])
        else:
            print('nothing found in the qualfile - maybe another GAMMA version?')
            print('please check this file manually: {}'.format(qualfile))
            return np.nan, np.nan, np.nan
    return daz_icc, dr_icc, daz_sd


def get_rangeshift_ICC(offile):
    # ok but this is 0 in offset files..
    rshift = float(grep1line('range_offset_polynomial', offile).split()[1])
    return rshift


def fix_oldorb_update_lt(ltfile, offile, azshiftm = 0.039):
    #samples = int(grep1line('interferogram_width', offile).split()[-1])
    #azoff = float(grep1line('azimuth_offset_polynomial', offile).split()[1])
    #offile2 = 
    # to correct:
    #gc_map_fine ltfile samples offile2 ltfile_ok 1
    
    ltfile_out = ltfile #+'.ok'
    lines = int(grep1line('interferogram_azimuth_lines', offile).split()[-1])
    samples = int(grep1line('interferogram_width', offile).split()[-1])
    azires = float(grep1line('interferogram_azimuth_pixel_spacing', offile).split()[1])
    a = np.fromfile(ltfile,dtype=np.complex64).byteswap()
    b = a.reshape((lines,samples))
    az = np.imag(b)
    rg = np.real(b)
    azshiftpx = azshiftm/azires
    # correct only non-zero values
    az[az!=0]-=azshiftpx
    # return to cpx
    cpx = rg + 1j*az
    cpx.byteswap().tofile(ltfile_out)


def fix_oldorb_update_off(offile, azshiftm=-0.039, returnval = False):
    '''
    offset file azimuth shift in azimuth_offset_polynomial appears to be in SLC pixel, not multilooked..
    ok - it is actually all correct (12/2022)
    '''
    azires = float(grep1line('interferogram_azimuth_pixel_spacing', offile).split()[1])
    azilooks = int(grep1line('interferogram_azimuth_looks', offile).split()[1])
    aziorigres = azires / azilooks
    # get previous estimate from off file
    azioffset = float(grep1line('azimuth_offset_polynomial', offile).split()[1])
    azshiftpx = azshiftm/aziorigres
    # should be -39 mm here to fit the RDC-resampled data (although directly using POD would give better estimate)
    aziok = azioffset + azshiftpx
    # get it back to the off file
    oldstr='azimuth_offset_polynomial.*'
    newstr='azimuth_offset_polynomial: '+str(aziok)+' 0.0 0.0 0.0 0.0 0.0'
    sed_replace(oldstr, newstr, offile)
    if returnval:
        return aziok



def fix_oldorb_shift_oneoff_track(track=2):
    ''' Careful - to be run only once!
    This will modify LUT tables and coreg_quality files, so that daz will include shift due to updated orbits.
    This means: for all O and OR, shift the values by -39 mm, where O=epoch resampled using old orbits, OR=epoch that was using O as RSLC3
    '''
    frames=os.listdir(os.path.join(os.environ['LiCSAR_procdir'], str(track)))
    for frame in frames:
        if len(frame) != 17:
            continue
        try:
            master = fc.get_master(frame)
        except:
            continue
        if not os.path.exists(os.path.join(os.environ['LiCSAR_public'], str(track), frame)):
            print('frame '+frame+' not properly initialized, deleting')
            os.system('rm -rf {0}'.format(os.path.join(os.environ['LiCSAR_procdir'], str(track), frame)))
            continue
        #masterfile = os.path.join(os.environ['LiCSAR_procdir'], track, frame, 'SLC', master, )
        #cdate_master = 
        if int(master) > 20210614:
            # this was day of updating orbits. everything above is OK!
            print('frame '+frame+' is ok')
            continue
        print('processing frame '+frame)
        fix_oldorb_shift_oneoff(frame)

'''
        eofile = glob.glob(os.path.join(os.environ['LiCSAR_procdir'], track, frame, 'SLC', master, '*EOF'))[0]
        eofile=os.path.basename(eofile)
        eofdate = eofile.split('_')[5].split('T')[0]
'''

import datetime as dt
import warnings
warnings.filterwarnings("ignore")
from LiCSquery import *
import shutil

def fix_oldorb_shift_oneoff(frame, tmpdir = '/work/scratch-pw3/licsar/earmla/temp3/'):
    ''' Careful - to be run only once!
    This will modify LUT tables and coreg_quality files, so that daz will include shift due to updated orbits.
    This means: for all O and OR, shift the values by -39 mm, where O=epoch resampled using old orbits, OR=epoch that was using O as RSLC3
    '''
    track = str(int(frame[:3]))
    framedir=os.path.join(os.environ['LiCSAR_procdir'], track, frame)
    bckdir=os.path.join(os.environ['LiCSAR_procdir'], track, frame, 'backup')
    tmpdir = tmpdir+frame
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    master = fc.get_master(frame)
    try:
        mastersat = fc.get_master(frame, asfilenames=True)[0][2]
    except:
        try:
            mastersat = os.path.basename(glob.glob(framedir+'/SLC/*/*EOF')[0])[2]
        except:
            mastersat = 'A'
    lutdir = os.path.join(os.environ['LiCSAR_procdir'], track, frame, 'LUT')
    logdir = os.path.join(os.environ['LiCSAR_procdir'], track, frame, 'log')
    #
    table = pd.DataFrame(columns=['epoch', 'mdate', 'RSLC3', 'azshift_SD', 'daz_SD', 'daz_ICC', 'dr_ICC'])
    print('checking and updating')
    if not os.path.exists(lutdir):
        os.mkdir(lutdir)
    luts = glob.glob(lutdir+'/20??????.7z')
    if len(luts)>0:
        for z in [ os.path.basename(a) for a in luts ]:
            epoch = z.split('.')[0]
            #if int(epoch) > 20200800:
            #print('processing epoch '+epoch)
            # first process the LUT file:
            lutfile = os.path.join(lutdir,z)
            mdate = dt.datetime.fromtimestamp(os.path.getmtime(lutfile))
            mdate = int(mdate.strftime('%Y%m%d'))
            #rc = os.system('cd {0}; 7za x {1} {2}/{3}_{2}.off>/dev/null'.format(tmpdir, lutfile, epoch, master))
            #ltfile = os.path.join(tmpdir, epoch, master+'_'+epoch+'.slc.mli.lt')
            offile = os.path.join(tmpdir, epoch, master+'_'+epoch+'.off')
            if not os.path.exists(offile):
                rc = os.system('cd {0}; 7za x {1} {2}/{3}_{2}.off>/dev/null'.format(tmpdir, lutfile, epoch, master))
            if os.path.exists(offile):
                try:
                    azishift_SD = get_azshift_SD(offile)
                except:
                    print('bad off file - deleting '+lutfile)
                    rc = os.system('rm '+lutfile)
                    azishift_SD = np.nan
            else:
                # error with LUT file, mv to bck:
                if not os.path.exists(bckdir):
                    os.mkdir(bckdir)
                rc=shutil.move(os.path.join(lutdir,epoch+'.7z'), os.path.join(bckdir,epoch+'.7z'))
                azishift_SD = np.nan
            # now get more info from qualfile
            fixrslc3 = False
            qualfile = os.path.join(logdir, 'coreg_quality_{0}_{1}.log'.format(master, epoch))
            if os.path.exists(qualfile):
                try:
                    daz_icc, dr_icc, daz_sd = get_shifts_from_qualfile(qualfile)
                except:
                    print('Some error reading qualfile - please check this file manually: {}'.format(qualfile))
                    daz_icc, dr_icc, daz_sd = np.nan, np.nan, np.nan
                    daz_sd = azishift_SD
                try:
                    rslc3 = grep1line('Spectral diversity estimation',qualfile).split(':')[-1]
                    if not rslc3:
                        fixrslc3 = True
                    else:
                        rslc3=int(rslc3.strip())
                except:
                    print('error in qualfile to extract RSLC3 info')
                    fixrslc3 = True
            else:
                print('ERROR - coreg qual file does not exist')
                daz_icc, dr_icc, daz_sd = np.nan, np.nan, np.nan
                daz_sd = azishift_SD
                fixrslc3 = True
            if fixrslc3:
                if np.abs((pd.Timestamp(master)-pd.Timestamp(epoch)).days)<180:
                    rslc3 = int(master)
                else:
                    try:
                        rslc3= int(mdate)
                    except:
                        rslc3= int(master)
            if azishift_SD == np.nan:
                azishift_SD = daz_sd
            if azishift_SD != np.nan:
                # table = table.append({'epoch':int(epoch), 'mdate': mdate, 'RSLC3': rslc3, 'azshift_SD': azishift_SD, 'daz_SD': daz_sd, 'daz_ICC': daz_icc, 'dr_ICC': dr_icc}, ignore_index=True)
                newpdline = pd.DataFrame({'epoch': [int(epoch)], 'mdate': [mdate], 'RSLC3': [rslc3], 'azshift_SD': [azishift_SD], 'daz_SD': [daz_sd],
                     'daz_ICC': [daz_icc], 'dr_ICC': [dr_icc]})
                table = pd.concat([table, newpdline], ignore_index=True)
    #
    table['epochdate'] = table.epoch.astype(int).astype(str)
    if not table.empty:
        table['epochdate'] = table.apply(lambda x : pd.to_datetime(str(x.epochdate)).date(), axis=1)
    dazdb = get_daz_frame(frame)
    # fill non-existing
    for i,row in table[table.isna().sum(axis=1)==0].iterrows():
        if not np.isin(row.epochdate, dazdb.epoch.values):
        #if row.epochdate not in dazdb.epoch:
            print('updating database for epoch '+str(row.epoch))
            epoch=str(row.epoch)
            rslc3=str(row.RSLC3)
            daz=row.azshift_SD
            ccazi=row.daz_ICC
            ccrg = row.dr_ICC
            try:
                eofile = glob.glob(os.path.join(os.environ['LiCSAR_procdir'], track, frame, 'SLC', master, '*EOF'))[0]
                orb=os.path.basename(eofile)
            except:
                orb='imported from LUT'
            lq.ingest_esd(frame, epoch, rslc3, daz, ccazi, ccrg, orb, overwrite=False)
    # add from non-fixed db records
    dazdb = dazdb[dazdb['orbfile']!='fixed_as_in_GRL']
    for i,row in dazdb.iterrows():
        if not np.isin(row.epoch, table.epochdate.values):
            print('adding from db :'+str(row.epoch))
            epoch = int(row.epoch.strftime('%Y%m%d'))
            try:
                mdate = int(row.orbfile.split('_')[5].split('T')[0])
            except:
                #mdate = np.nan
                if epoch > int(master):
                    mdate = epoch   # assumming here only - perhaps better than with nan
                else:
                    mdate = int(master)
            rslc3 = int(row.rslc3.strftime('%Y%m%d'))
            azishift_SD = row.daz
            daz_sd = row.daz
            daz_icc = row.cc_azi
            dr_icc = row.cc_range
            newpdline = pd.DataFrame({'epoch':[int(epoch)], 'mdate': [mdate],
            'RSLC3': [rslc3], 'azshift_SD': [azishift_SD], 'daz_SD': [daz_sd],
            'daz_ICC': [daz_icc], 'dr_ICC': [dr_icc], 'epochdate': [row['epoch']]})
            table = pd.concat([table, newpdline], ignore_index=True)
        # get epochs 'to correct', i.e. that used old orb files, and if they were used for RSLC3:
        #tocorrectepochs = []
    O = table[table['epoch']<20200729]
    O = O[O['mdate']<20210614]
    tocorrectepochs = O.epoch.astype(int).astype(str)
    if not tocorrectepochs.empty:
        tocorrectepochs = list(tocorrectepochs.values)
        checkit = 1
    else:
        checkit = 0
    if checkit == 1:
        allepochs = list(table.epoch.astype(int).values) + [int(master)]
        missingrslc3s = table[~table.RSLC3.astype(int).isin(allepochs)].RSLC3.astype(int).unique()
        if len(missingrslc3s)>0:
            print('missing rslcs detected, substituting')
            #possible_rslc3s = table[table.RSLC3.astype(int).isin(allepochs)].RSLC3.astype(int).unique()
            affected = table[table.RSLC3.astype(int).isin(missingrslc3s)].epoch.astype(int).values
            #possible_rslc3s = table[~table.epoch.astype(int).isin(affected)].epoch.astype(int).values
            possible_rslc3s = table[~table.epoch.astype(int).isin(affected)].epochdate #.values
            if not possible_rslc3s.empty:
                possible_rslc3s_ord = possible_rslc3s.apply(lambda x : pd.to_datetime(str(x)).toordinal()).values
                possible_rslc3s = table[~table.epoch.astype(int).isin(affected)].epoch.values
                for missing in missingrslc3s:
                    agg=pd.Timestamp(str(missing)).toordinal()
                     #.unique()
                    #possible_rslc3s = table[table.RSLC3.astype(int).isin(allepochs)].RSLC3.astype(int).unique()
                    substitute_i = np.argmin(np.abs(possible_rslc3s_ord-agg))
                    substitute = possible_rslc3s[substitute_i]
                    selec = table[table['RSLC3'].astype(int)==missing].index
                    table.loc[selec,'RSLC3'] = substitute
            #
    while checkit == 1:
        tocheck = table[table['RSLC3'].astype(int).astype(str).isin(tocorrectepochs)]
        #tocheck = tocheck[~tocheck['epoch'].astype(int).astype(str).isin(tocorrectepochs)]
        tocheck = tocheck[~tocheck['epoch'].astype(int).astype(str).isin(tocorrectepochs)].epoch.astype(int).astype(str)
        if tocheck.empty:
            checkit = 0
        else:
            print('iteration of OR')
            tocorrectepochs = tocorrectepochs + list(tocheck.values)
            #
    if len(tocorrectepochs)>0:
        tocorrectepochs = list(set(tocorrectepochs)) # remove duplicates
        tocorrectepochs.sort()
        print('correcting '+str(len(tocorrectepochs))+' epochs between '+str(tocorrectepochs[0])+' and '+str(tocorrectepochs[-1]))
        if not os.path.exists(bckdir):
            os.mkdir(bckdir)
    for epoch in tocorrectepochs:
        # now for those needed, update the value in off:
        offile = os.path.join(tmpdir, epoch, master+'_'+epoch+'.off')
        if os.path.exists(os.path.join(lutdir,epoch+'.7z')):
            try:
                rc=shutil.copyfile(os.path.join(lutdir,epoch+'.7z'), os.path.join(bckdir,epoch+'.7z'))
            except:
                print('error copying LUT of '+str(epoch))
        try:
            try:
                newazishift = fix_oldorb_update_off(offile, azshiftm=-0.039, returnval = True)
                rc = os.system('cd {0}; 7za u {1} {2}/{3}_{2}.off>/dev/null'.format(tmpdir, os.path.join(lutdir,epoch+'.7z'), epoch, master))
            except:
                print('error updating off file in LUT of '+str(epoch))
                newazishift = float(table[table.epoch == int(epoch)].azshift_SD)-39/14000
            # now also change the coreg_qual file - or just .. move it away...
            qualfile = os.path.join(logdir, 'coreg_quality_{0}_{1}.log'.format(master, epoch))
            rc = os.system('mv {0} {1}/.'.format(qualfile, bckdir))
            # and finally update in database!
            rc = update_esd(frame, epoch, colupdate = 'daz', valupdate = newazishift)
        except:
            print('some error with frame '+frame+', epoch '+str(epoch))


def fix_oldorb_pds(framespd, esds):
    """ Gets azimuth correction for all frames in framespd. Unfinished..
    """
    for frame in framespd.frame.values:
        dazes = get_daz_frame(frame)
        epochs = []
        epochs = epochs + dazes[dazes['orbfile']==''].epoch.to_list()
        epochs = epochs + dazes[dazes['orbfile']=='fixed_as_in_GRL'].epoch.to_list()  # this needs to be recorrected (+39)
        # finally, need to find epochs where the correction is as rslc3, and correct them as well, plus other in the cascade.....
        if epochs:
            azispd = get_azioffs_old_new_POD(frame, epochs)
    return esds


def flag_old_new_POD(esds):
    """ This will add new column (boolean): 'new_POD' where True means 'safe to use' as it uses orbits >2020/07/30"""
    esds['new_POD'] = False
    for frame, group in esds.groupby('frame'):
        # first check if there is any epoch to fix (maybe not?)
        dazes = get_daz_frame(frame)
        dazes = dazes[np.isin(dazes['epoch'], group['epochdate'])]
        epochs = []
        #epochs = epochs + dazes[dazes['orbfile']==''].epoch.to_list()
        epochs = epochs + dazes[dazes['orbfile']=='fixed_as_in_GRL'].epoch.to_list()
        aa = pd.Series(epochs)
        epochs = aa[aa < dt.datetime(2020, 7, 30).date()].to_list() # sometimes we have no info on POD used in the db..
        lenep = 0
        while len(epochs)>lenep:
            epochs = epochs + dazes[np.isin(dazes['rslc3'], epochs)].epoch.to_list()
            epochs = list(set(epochs))
            lenep = len(epochs)
        group.loc[~np.isin(group['epochdate'], epochs), 'new_POD'] = True
        esds.update(group)
    return esds


def get_azioffs_old_new_POD(frame, epochs = None):
    """ Function to get correction for PODs established after end of July 2020 in azimuth
    """
    print('getting old/new POD difference corrections for frame '+frame)
    datelim = dt.datetime(2020,7,31).date()
    if type(epochs) == type(None):
        epochs = fc.get_epochs(frame, return_as_dt=True) #2018-09-01
    master_s1ab = get_frame_master_s1ab(frame)
    master = fc.get_master(frame, asdatetime = True)
    azioffs = []
    selepochs = []
    for epoch in epochs:
        if epoch > datelim:
            print('epoch ' + str(epoch) + ' was surely processed with POD v1.4+, skipping')
            continue
        epoch_s1ab = flag_s1b([epoch], master, master_s1ab, returnstr=True )[0]
        timesample = dt.datetime.combine(epoch, master.time())
        neworbs = get_orbit_filenames_for_datetime(timesample, 'POEORB', s1ab='S1'+epoch_s1ab)
        oldorbs = getoldorbpath(neworbs)
        neworb = neworbs[0]
        oldorb = oldorbs[0]
        if not oldorb:
            print('none old for '+str(epoch))
            continue
        azioff = get_azi_diff_from_two_orbits(oldorb, neworb, timesample)
        selepochs.append(epoch)
        azioffs.append(azioff)
    if not selepochs:
        print('no epoch was selected for correction')
        return False
    selazis = np.array(azioffs) *1000
    azispd = pd.DataFrame({'epochdate': selepochs,
     'pod_diff_azi_mm': selazis})
    return azispd



def get_azshift_lt(ltfile = '20210425.mli.lt', offile = '20210413_20210425.off.start', az_ml = 4, rg_ml = 20, return_rg = True):
    lines = int(grep1line('interferogram_azimuth_lines', offile).split()[-1])
    samples = int(grep1line('interferogram_width', offile).split()[-1])
    #azi_res = float(grep1line('interferogram_azimuth_pixel_spacing', offile).split()[1])
    #azi_full_lines = az_ml*lines
    
    a = np.fromfile(ltfile,dtype=np.complex64).byteswap()
    b = a.reshape((lines,samples))
    az = np.imag(b)
    # get az shift at the middle range:
    midsample = int(np.floor(samples/2))
    m = az[midsample,:]
    m = m[m != 0]
    azshift = az_ml * (np.mean(m) - midsample)
    #azshift_SD = float(grep1line('azimuth_offset_polynomial', offile).split()[1])
    #azshift = azshift + azshift_SD
    if not return_rg:
        return azshift
    else:
        rg = np.real(b)
        midline = int(np.floor(lines/2))
        m = rg[:,midline]
        m = m[m != 0]
        rgshift = rg_ml * (np.mean(m) - midline)
        return azshift, rgshift


def get_gacos_in_coord(lon,lat,epochstr, frame, inmm=True, domean=True):
    #import rioxarray
    epochpath = os.path.join(os.environ['LiCSAR_public'], str(int(frame[:3])), frame, 'epochs', epochstr)
    gacospath = os.path.join(epochpath, epochstr+'.sltd.geo.tif')
    if not os.path.exists(gacospath):
        return np.nan
    else:
        try:
            f = rioxarray.open_rasterio(gacospath)
            if domean:
                out=float(f[0].mean())
            else:
                out = float(f[0].sel(x=lon, y=lat, method='nearest').values)
        except:
            print('error reading')
            print(gacospath)
            return np.nan
        if inmm:
            out = rad2mm_s1(out)
        return out


def get_table_azishifts(frame):
    track = str(int(frame[:3]))
    #frame = '099A_05416_131313'
    #frame = '099A_05417_131313'
    master = fc.get_master(frame)
    tmpdir = '/work/scratch-pw3/licsar/earmla/temp'
    lutdir = os.path.join(os.environ['LiCSAR_procdir'], track, frame, 'LUT')
    logdir = os.path.join(os.environ['LiCSAR_procdir'], track, frame, 'log')
    table = pd.DataFrame(columns=['epoch', 'azshift_RDC_ICC', 'rgshift_RDC_ICC', 'azshift_SD', 'daz_ICC', 'dr_ICC'])
    rc = os.system('rm -rf {0}/*'.format(tmpdir))
    for z in os.listdir(lutdir):
        epoch = z.split('.')[0]
        #if int(epoch) > 20200800:
        print('processing epoch '+epoch)
        qualfile = os.path.join(logdir, 'coreg_quality_{0}_{1}.log'.format(master, epoch))
        if os.path.exists(qualfile):
            try:
                daz_icc, dr_icc, daz_sd = get_shifts_from_qualfile(qualfile)
            except:
                print('Some error reading qualfile - please check this file manually: {}'.format(qualfile))
                daz_icc, dr_icc, daz_sd = np.nan, np.nan, np.nan
        else:
            print('ERROR - coreg qual file does not exist')
            daz_icc, dr_icc, daz_sd = np.nan, np.nan, np.nan
        rc = os.system('cd {0}; 7za x {1} >/dev/null'.format(tmpdir, os.path.join(lutdir,z)))
        ltfile = os.path.join(tmpdir, epoch, master+'_'+epoch+'.slc.mli.lt')
        offile = os.path.join(tmpdir, epoch, master+'_'+epoch+'.off')
        if os.path.exists(ltfile) and os.path.exists(offile):
            azishift, rgshift = get_azshift_lt(ltfile, offile)
            azishift_SD = get_azshift_SD(offile)
            newpdline = pd.DataFrame({'epoch':[epoch], 'azshift_RDC_ICC': [azishift], 'rgshift_RDC_ICC': [rgshift],
                                      'azshift_SD': [azishift_SD], 'daz_ICC': [daz_icc], 'dr_ICC': [dr_icc]})
            table = pd.concat([table, newpdline], ignore_index=True)
        else:
            print('files not extracted correctly, skipping epoch '+epoch)
        rc = os.system('rm -rf {0}/*'.format(tmpdir))
    table['epochdate'] = table.apply(lambda x : pd.to_datetime(str(x.epoch)).date(), axis=1)
    return table



def test_rdcs():
    '''
    this shows effect of RDC_TRANS on the azimuth shifts
    '''
    frame = '099A_05416_131313'
    table = get_table_azishifts(frame)
    #or: 
    #table = pd.read_csv('099A_new_with_all.csv')
    table['epochdate'] = table.apply(lambda x : pd.to_datetime(str(x.epoch)).date(), axis=1)
    
    table_old = table.set_index('epochdate')
    frame = '099A_05417_131313'
    table = get_table_azishifts(frame)
    table_new = table.set_index('epochdate')
    
    az_resolution = 13.9751 # m
    cols = ['azshift_RDC_ICC', 'azshift_SD', 'daz_ICC'] #, 'daz_SD']
    oldshifts = table_old[cols]*az_resolution*1000
    newshifts = table_new[cols]*az_resolution*1000
    
    rg_resolution = 2.329562 # m
    col = 'dr_ICC'
    oldshifts[col] = table_old[col]*rg_resolution*1000
    newshifts[col] = table_new[col]*rg_resolution*1000
    col = 'rgshift_RDC_ICC'
    oldshifts[col] = table_old[col]*rg_resolution*1000
    newshifts[col] = table_new[col]*rg_resolution*1000
    
    #oldshifts = oldshifts
    newshifts.to_csv('099A_new_incl_rg.csv')
    oldshifts.to_csv('099A_old_incl_rg.csv')
    #table = pd.read_csv('099A_new_with_all.csv')
    table['epochdate'] = table.apply(lambda x : pd.to_datetime(str(x.epochdate)).date(), axis=1)
    
    table_new['azshift_RDC'] = table_new['azshift_RDC_ICC'] - table_new['daz_ICC']
    table_old['azshift_RDC'] = table_old['azshift_RDC_ICC'] - table_old['daz_ICC']
    table_new['azshift_RDC_ICC_SD'] = table_new['azshift_RDC_ICC'] + table_new['azshift_SD']
    table_old['azshift_RDC_ICC_SD'] = table_old['azshift_RDC_ICC'] + table_old['azshift_SD']
    
    difft = table_new - table_old
    difft = difft[~np.isnan(difft['daz_ICC_SD'])]
    difft = difft[difft['azshift_RDC_ICC']>-200]
    difft = difft[difft['azshift_RDC_ICC']<200]
    difft = difft[difft['azshift_SD']>-200]
    difft = difft[difft['azshift_SD']<200]
    
    difft['daz_ICC_SD_noTI'] = esds_new.daz_mm_notide_noiono_grad - esds_old.daz_mm_notide_noiono_grad
    
    ax1 = difft['azshift_RDC_ICC_SD'].plot()
    ax1.set_ylim(-150,50)
    plt.show()
    
    oldshifts.plot()
    newshifts.plot()
    plt.show()
    
    diff = newshifts - oldshifts
    diff = diff[diff<80]
    diff.plot()
    plt.show()

