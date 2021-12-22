
from daz_lib import *
import os

#######################################
# step 1 - export ESDs to a txt file:
# see prepare_daz_licsar.sh
# or just get it from licsinfo's esd table..
# (see LiCSquery - get_daz etc.)

dazfile = 'esds.txt'

#######################################
# step 2 - add additional information to the esds.txt and frames.txt files
# to be run at JASMIN
# this is to prepare the framespd and esds

fname = 'esds2021_frames.txt'
outcsv = 'framespd_2021.csv'
if not os.path.exists(outcsv):
    generate_framespd(fname, outcsv)


#######################################
# step 3 - get solid Earth tides
################### SOLID EARTH TIDES
# see get_SET.sh


#######################################
# step 4 - merge and process further in python
###########################
# now (finally) load esds and tides to python, merge etc.
framescsv = 'framespd_2021.csv'
esdscsv = 'esds2021_take4.txt'
tidescsv = 'earthtides_take4.csv'

earthtides = pd.read_csv(tidescsv)
esds = pd.read_csv(esdscsv)
framespd = pd.read_csv(framescsv)

esds = merge_tides(esds, framespd, earthtides)

esds['epochdate'] = esds.apply(lambda x : pd.to_datetime(str(x.epochdate)).date(), axis=1)
esds = esds.drop('epoch', axis=1)
mindate = esds['epochdate'].min()
maxdate = esds['epochdate'].max()


#prepare esds:
esds, framespd = df_preprepare_esds(esds, framespd, mindate)
esds = esds.reset_index(drop=True)
framespd = framespd.reset_index(drop=True)

# to back up before continuing:
#framespd.to_csv('backup_0623_framespd.csv')
#esds.to_csv('backup_0623_esds.csv')


#######################################
# step 5 - get daz iono
################### IONOSPHERE 


# estimating the ionosphere - takes long (several hours)
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
        daz_iono_grad, hionos, tecs_A_master, tecs_B_master = calculate_daz_iono(frame, esds, framespd, method = 'gomba', out_hionos = True, out_tec_master = True)
        hiono = np.mean(hionos)
        hiono_std = np.std(hionos)
    except:
        print('some error occurred here')
        continue
    esds.at[esds[esds['frame']==frame].index, 'daz_iono_grad_mm'] = daz_iono_grad*resolution*1000
    #esds.at[esds[esds['frame']==frame].index, 'daz_iono_with_F2'] = daz_iono_with_F2
    esds.at[esds[esds['frame']==frame].index, 'daz_mm_notide_noiono_grad'] = esds[esds['frame']==frame]['daz_mm_notide'] - esds[esds['frame']==frame]['daz_iono_grad_mm'] #*resolution*1000
    framespd.at[framespd[framespd['frame']==frame].index, 'Hiono'] = hiono
    framespd.at[framespd[framespd['frame']==frame].index, 'Hiono_std'] = hiono_std
    framespd.at[framespd[framespd['frame']==frame].index, 'Hiono_range'] = max(hionos)-min(hionos)
    framespd.at[framespd[framespd['frame']==frame].index, 'tecs_A'] = tecs_A_master
    framespd.at[framespd[framespd['frame']==frame].index, 'tecs_B'] = tecs_B_master
    #esds.at[esds[esds['frame']==frame].index, 'daz_mm_notide_noiono_F2'] = esds[esds['frame']==frame]['daz_mm_notide'] - esds['daz_iono_with_F2']*resolution*1000


# to back up before continuing:
framespd.to_csv('backup_0921_framespd_iono_check.csv')
esds.to_csv('backup_0921_esds_iono_check.csv')


#######################################
# step 6 - get plate motion model
################### ITRF2014


#takes again long - some 2-3 hours
print('getting ITRF2014 values - using the average for the 222x222 km around the frame centre')
framespd = df_get_itrf_slopes(framespd)

framespd.to_csv('backup_1907_framespd_iono_itrf.csv')


# step 6+ -- correct for daz_ARP=-39 mm
#################
cols = ['daz_mm','daz_mm_notide', 'daz_mm_notide_noiono_grad']
ep = esds[esds.epochdate > pd.Timestamp('2020-07-30')][cols]
esds.update(ep.subtract(-39))





# 2021-10-12: the original way, only this time with full dataset:
esds, framespd = df_calculate_slopes(esds, framespd, alpha = 1, eps = 1.35, bycol = 'daz_mm', subset = False, roll_assist = True)
esds, framespd = df_calculate_slopes(esds, framespd, alpha = 1, eps = 1.35, bycol = 'daz_mm_notide', subset = False, roll_assist = True)
esds, framespd = df_calculate_slopes(esds, framespd, alpha = 1, eps = 1.35, bycol = 'daz_mm_notide_noiono_grad', subset = False, roll_assist = True)
#esds, framespd = df_calculate_slopes(esds, framespd, alpha = 1, eps = 1.35, bycol = 'daz_mm_notide_noiono_F2')

# to back up before continuing:
framespd.to_csv('backup_1012_framespd_iono_full_with_slopes_nosubset_roll.csv')
esds.to_csv('backup_1012_esds_iono_full_with_slopes_nosubset_roll.csv')

esds_roll = esds.copy(deep=True)
framespd_roll = framespd.copy(deep=True)

# skipping roll assistance - as e.g. Turkey is not correct in ITR2014 PMM
esds, framespd = df_calculate_slopes(esds, framespd, alpha = 1, eps = 1.35, bycol = 'daz_mm', subset = False, roll_assist = False)
esds, framespd = df_calculate_slopes(esds, framespd, alpha = 1, eps = 1.35, bycol = 'daz_mm_notide', subset = False, roll_assist = False)
esds, framespd = df_calculate_slopes(esds, framespd, alpha = 1, eps = 1.35, bycol = 'daz_mm_notide_noiono_grad', subset = False,roll_assist = False)
#esds, framespd = df_calculate_slopes(esds, framespd, alpha = 1, eps = 1.35, bycol = 'daz_mm_notide_noiono_F2')

# to back up before continuing:
framespd.to_csv('backup_1012_framespd_iono_full_with_slopes_nosubset_noroll.csv')
esds.to_csv('backup_1012_esds_iono_full_with_slopes_nosubset_noroll.csv')

