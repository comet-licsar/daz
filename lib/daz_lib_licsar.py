#!/usr/bin/env python3

import LiCSquery as lq
import numpy as np
from LiCSAR_lib.LiCSAR_misc import *
import os, glob
import pandas as pd
import framecare as fc
try:
    from orbit_lib import *
except:
    print('LiCSAR orbit_lib not found, cannot process orbit files')


from daz_lib import *


# e.g.
# framelist=pd.read_csv('frames.txt'); framelist=list(framelist.frame)
# where you may get frames.txt e.g. by:
# cd $LiCSAR_procdir
# for tr in `seq 1 175`; do for f in `ls $tr`; do m=`ls $tr/$f/SLC | head -n1`; hgtfile=$LiCSAR_public/$tr/$f/metadata/$f'.geo.hgt.tif'; ll=`gdalinfo $hgtfile | grep ^Center`; lon=`echo $ll | cut -d "," -f1 | cut -d '(' -f2`; lat=`echo $ll | cut -d ")" -f1 | cut -d ',' -f2`; echo $f","$m","$lon","$lat >> $outfr;  done;done
def extract2txt_esds_all_frames(framelist, outfile='esds.txt'):
    dazes=pd.DataFrame()
    for frame in framelist:
        try:
            a=extract2txt_esds_frame(frame)
            dazes=dazes.append(a)
        except:
            print('frame '+frame+' is empty')
    dazes=dazes.reset_index(drop=True)
    dazes.to_csv(outfile, index=False)


def extract2txt_esds_frame(frame):
    '''
    extracts to esds txt full data for given frame, from database
    the resultant txt file is a csv as:
    frame,esd_master,epoch,daz_total_wrt_orbits,daz_cc_wrt_orbits,orbits_precision,version
    '''
    a = get_daz_frame(frame)
    a['epoch']=a.epoch.apply(lambda x: x.strftime('%Y%m%d'))
    a['esd_master']=a.rslc3.apply(lambda x: x.strftime('%Y%m%d'))
    a['daz_total_wrt_orbits']=a.daz+a.cc_azi
    a['orbits_precision'] = 'P'  # only Ps should be in database
    a['version'] = 'm' # i forgot what this is for, but should be ok any letter (?)
    a['frame'] = frame
    a=a.rename(columns={'cc_azi':'daz_cc_wrt_orbits'})
    return a[['frame','esd_master','epoch','daz_total_wrt_orbits','daz_cc_wrt_orbits','orbits_precision','version']]


def get_daz_frame(frame):
    polyid=lq.get_frame_polyid(frame)[0][0]
    daztb = lq.do_pd_query('select * from esd where polyid={};'.format(polyid))
    return daztb


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


def get_frame_master_s1ab(frame):
    tr = int(frame[:3])
    metafile = os.path.join(os.environ['LiCSAR_public'], str(tr), frame, 'metadata', 'metadata.txt')
    if not os.path.exists(metafile):
        print('metadata file does not exist for frame '+frame)
        return 'X'
    primepoch = grep1line('master=',metafile).split('=')[1]
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
    #a['kr']=0.00
    a['dfDC'] = 0.00
    a['avg_height'] = 0.00
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
    a = extract_frame_master_s1abs(a)
    a.to_csv(outcsv, float_format='%.4f', index=False)



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
    if not returnperswath:
        dfDC = np.mean(dfDC)
        ka = np.mean(kas)
    #kr = np.mean(krs)
    if returnka:
        return dfDC, ka #, kr
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
    # should be -39 mm here to fit the RDC-resampled data.. checking now
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
                table = table.append({'epoch':int(epoch), 'mdate': mdate, 'RSLC3': rslc3, 'azshift_SD': azishift_SD, 'daz_SD': daz_sd, 'daz_ICC': daz_icc, 'dr_ICC': dr_icc}, ignore_index=True)
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
            table = table.append({'epoch':int(epoch), 'mdate': mdate, 
            'RSLC3': rslc3, 'azshift_SD': azishift_SD, 'daz_SD': daz_sd, 
            'daz_ICC': daz_icc, 'dr_ICC': dr_icc, 'epochdate': row.epoch}, ignore_index=True)
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


def get_azioffs_old_new_POD(frame, epochs = None):
    """ Function to get correction for PODs established after 2020-07-31 in azimuth
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
            break
        epoch_s1ab = flag_s1b([epoch], master, master_s1ab, returnstr=True )[0]
        timesample = dt.datetime.combine(epoch, master.time())
        neworbs = get_orbit_filenames_for_datetime(timesample, 'POEORB', s1ab='S1'+epoch_s1ab)
        oldorbs = getoldorbpath(neworbs)
        neworb = neworbs[0]
        oldorb = oldorbs[0]
        if not oldorb:
            print('none old for '+str(epoch))
            break
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
            table = table.append({'epoch':epoch, 'azshift_RDC_ICC': azishift, 'rgshift_RDC_ICC': rgshift, 'azshift_SD': azishift_SD, 'daz_ICC': daz_icc, 'dr_ICC': dr_icc}, ignore_index=True)
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

