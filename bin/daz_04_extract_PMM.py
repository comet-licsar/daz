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

Outputs :
 - frames_with_itrf.csv - new columns with the ITRF2014 PMM in azimuth direction

=====
Usage
=====
daz_04_extract_PMM.py [--infra frames.csv] [--outfra frames_with_itrf.csv]

"""
#%% Change log
'''
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
    
    #%% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help", "infra", "outfra"])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == "--infra":
                inframesfile = a
            elif o == "--outfra":
                outframesfile = a
        
        if os.path.exists(outframesfile):
            raise Usage('output frames csv file already exists. Cancelling')
        if not os.path.exists(inframesfile):
            raise Usage('input frames csv file does not exist. Cancelling'))
            
    except Usage as err:
        print("\nERROR:",)
        print("  "+str(err.msg))
        print("\nFor help, use -h or --help.\n")
        return 2
    
    framespd = pd.read_csv(inframesfile)
    
    # get plate motion model
    ################### ITRF2014
    #takes again long - some 2-3 hours
    print('getting ITRF2014 values - using the average for the 222x222 km around the frame centre')
    framespd = df_get_itrf_slopes(framespd)
    framespd.to_csv(outframesfile)
    print('done')

#%% main
if __name__ == "__main__":
    sys.exit(main())



