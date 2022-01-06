#!/usr/bin/env python3
"""
v1.0 2022-01-03 Milan Lazecky, Leeds Uni

This script will decompose azimuth velocities of frames into N,E component in a regular global grid.
It includes RMSE_VEL estimation as a weighted LS error term.

===============
Input & output files
===============
Inputs :
 - frames_final.csv - must contain results of velocity estimates (daz_05)

Outputs :
 - decomposed.csv

=====
Usage
=====
daz_06_decompose.py [--infra frames_final.csv] [--outdec decomposed.csv]

"""
#%% Change log
'''
v1.0 2022-01-03 Milan Lazecky, Uni of Leeds
 - Original implementation - based on codes from 2021-06-24
'''
# decompose
##########


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
    inframesfile = 'frames_final.csv'
    outdecfile = 'decomposed.csv'
    #%% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help", "infra =", "outdec ="])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == "--infra":
                inframesfile = a
            elif o == "--outdec":
                outdecfile = a
        
        if os.path.exists(outdecfile):
            raise Usage('output decomposition file already exists. Cancelling')
        if not os.path.exists(inframesfile):
            raise Usage('input frames csv file does not exist. Cancelling')
            
    except Usage as err:
        print("\nERROR:",)
        print("  "+str(err.msg))
        print("\nFor help, use -h or --help.\n")
        return 2
    
    # processing itself:
    framespd = pd.read_csv(inframesfile)
    print('decomposing frames')
    gridagg = decompose_framespd(framespd)
    print('getting ITRF 2014 PMM for new cells')
    gridagg = get_itrf_EN(gridagg)
    
    print('exporting final decomposed data to '+outdecfile)
    gridagg.to_csv(outdecfile)
    print('done')

#%% main
if __name__ == "__main__":
    sys.exit(main())



