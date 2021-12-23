#!/usr/bin/env python3

from daz_lib import *
from daz_timeseries import *


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

