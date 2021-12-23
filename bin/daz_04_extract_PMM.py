#!/usr/bin/env python3

from daz_lib import *

# get plate motion model
################### ITRF2014


#takes again long - some 2-3 hours
print('getting ITRF2014 values - using the average for the 222x222 km around the frame centre')
framespd = df_get_itrf_slopes(framespd)

framespd.to_csv('backup_1907_framespd_iono_itrf.csv')
