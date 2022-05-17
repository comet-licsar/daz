#!/usr/bin/env python3

import LiCSquery as lq
import numpy as np
from LiCSAR_lib.LiCSAR_misc import *
import os
import pandas as pd
import framecare as fc

# maybe not needed?
from daz_lib import *

def get_daz(frame):
    polyid=lq.get_frame_polyid(frame)[0][0]
    daztb = lq.do_pd_query('select * from esd where polyid={};'.format(polyid))


'''
how to get s1a/b? try this:
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


def fix_oldorb_update_lt(ltfile, offile, azshiftm = 0.0039):
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


def fix_oldorb_update_off(offile, azshiftm=-0.0039):
    '''
    offset file azimuth shift in azimuth_offset_polynomial appears to be in SLC pixel, not multilooked..
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
    tmpdir = '/work/scratch-pw/earmla/LiCSAR_temp/batchdir/temp2'
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

