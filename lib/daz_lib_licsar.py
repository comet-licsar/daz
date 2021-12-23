import LiCSquery as lq
import numpy as np
from LiCSAR_lib.LiCSAR_misc import *
import os
import pandas as pd
import framecare as fc

def get_daz(frame):
    polyid=lq.get_frame_polyid(frame)[0][0]
    daztb = lq.do_pd_query('select * from esd where polyid={};'.format(polyid))


def get_azshift_SD(offile):
    azshift_SD = float(grep1line('azimuth_offset_polynomial', offile).split()[1])
    return azshift_SD


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
    azshift_SD = float(grep1line('azimuth_offset_polynomial', offile).split()[1])
    return azshift_SD


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

