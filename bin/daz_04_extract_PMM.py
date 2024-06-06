#!/usr/bin/env python3
"""
v1.0 2022-01-03 Milan Lazecky, Leeds Uni

This script will extract plate motion model values from ITRF2014 PMM (extracted from UNAVCO website).
It will extract values in a buffer around given frame centre coordinate, and then average them.

===============
Input & output files
===============
Inputs :
 - frames.csv - either orig frames.csv or the one including iono columns
[- vel_gps.nc - GPS velocities, must contain VEL_E, VEL_N variables, expected in NNR]

Outputs :
 - frames_with_itrf.csv - new columns with the ITRF2014 PMM in azimuth direction

=====
Usage
=====
daz_04_extract_PMM.py [--add_eu] [--infra frames.csv] [--outfra frames_with_itrf.csv] [--velnc vel_gps_kreemer.nc]

Note: param velnc is optional, but if provided as nc file with VEL_E, VEL_N variables, it will be used as GPS velocities.
--add_eu will extract also ITRF2014 PMM EU along-track direction (ATD) that can be then used to transform the ATD velocities.
"""
#%% Change log
'''
v1.1 2024-06 ML
 - add EU in ATD
v1.0 2022-01-03 Milan Lazecky, Uni of Leeds
 - Original implementation - based on codes from 2021-06-24
'''

import getopt, os, sys
from daz_lib import *

class Usage(Exception):
    """Usage context manager"""
    def __init__(self, msg):
        self.msg = msg


#%% Main
def main(argv=None):
    
    #%% Check argv
    if argv == None:
        argv = sys.argv
    
    #%% Set default
    inframesfile = 'frames_with_iono.csv'
    outframesfile = 'frames_with_itrf.csv'
    velnc = 'vel_gps_kreemer.nc'
    add_eu = False

    #%% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help", "add_eu", "infra=", "outfra=", "velnc="])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == "--add_eu":
                add_eu = True
            elif o == "--infra":
                inframesfile = a
            elif o == "--outfra":
                outframesfile = a
            elif o == "--velnc":
                velnc = a
        
        if os.path.exists(outframesfile):
            raise Usage('output frames csv file already exists. Cancelling')
        if not os.path.exists(inframesfile):
            raise Usage('input frames csv file does not exist. Cancelling')
            
    except Usage as err:
        print("\nERROR:",)
        print("  "+str(err.msg))
        print("\nFor help, use -h or --help.\n")
        return 2
    
    framespd = pd.read_csv(inframesfile)
    
    # get plate motion model
    ################### ITRF2014
    #takes again long - some 2-3 hours
    print('getting plate motion model values using the average for the 222x222 km around the frame centre')
    print('(using ITRF2014 and external nc file for GPS, if available)')
    #framespd = df_get_itrf_slopes(framespd)
    framespd = df_get_itrf_gps_slopes(framespd, velnc=velnc, add_eu = add_eu)
    if add_eu:
        print('Finished extracting both NNR and EU PMM values. Note:')
        print("framespd['vel_eur'] = framespd['slope_from_daz'] - framespd['slope_vel_itrf_nnr'] + framespd['slope_vel_itrf_eu']")
    framespd.to_csv(outframesfile)
    print('done')

#%% main
if __name__ == "__main__":
    sys.exit(main())


'''
some notes that might be relevant:
gridagg['GPS_N_2014_EU'] = gridagg['GPS_N'] - gridagg['ITRF_N_2008'] + gridagg['ITRF_2014_EU_N']
gridagg['GPS_E_2014_EU'] = gridagg['GPS_E'] - gridagg['ITRF_E_2008'] + gridagg['ITRF_2014_EU_E']
framespd['vel_eur'] = framespd['slope_from_daz'] - framespd['slope_vel_itrf_nnr'] + framespd['slope_vel_itrf_eu']
'''
