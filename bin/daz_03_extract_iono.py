#!/usr/bin/env python3

from daz_lib import *
from daz_iono import *

# get daz iono
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
