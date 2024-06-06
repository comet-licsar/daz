#!/usr/bin/env python3
"""
v1.0 2022-01-03 Milan Lazecky, Leeds Uni

This script will generate kmz from frame time series

===============
Input & output files
===============
Inputs :
 - frames_final.csv - should have slopes calculated
 - esds_final.csv - should have iono corr.

Outputs :
 - esds.kmz

=====
Usage
=====
daz_export2kmz.py [--indaz esds_final.csv] [--infra frames_final.csv] [--outkmz esds.kmz]

"""
#%% Change log
'''
v1.0 2022-01-03 Milan Lazecky, Uni of Leeds
 - Original implementation - based on codes from 2021-06-24
######### update 2021-05-31: hopefully improved huber regression...

'''
from daz_lib import load_csvs
from daz_plotting import *

import getopt, os, sys

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
    indazfile = 'esds_final.csv'
    inframesfile = 'frames_final.csv'
    outkmzfile = 'esds.kmz'
    
    #%% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help", "indaz=", "infra=", "outkmz="])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == "--indaz":
                indazfile = a
            elif o == "--infra":
                inframesfile = a
            elif o == "--outkmz":
                outkmzfile = a
        
        if os.path.exists(outkmzfile):
            raise Usage('output esds kmz file already exists. Cancelling')
        if not os.path.exists(inframesfile):
            raise Usage('input frames csv file does not exist. Cancelling')
        if not os.path.exists(indazfile):
            raise Usage('input esds csv file does not exist. Cancelling')
            
    except Usage as err:
        print("\nERROR:",)
        print("  "+str(err.msg))
        print("\nFor help, use -h or --help.\n")
        return 2
    
    # processing itself:
    #esds = pd.read_csv(indazfile)
    #framespd = pd.read_csv(inframesfile)
    esds, framespd = load_csvs(esdscsv = indazfile, framescsv = inframesfile)
    if 'daz_mm_final' in esds:
        print('final values detected (2024/06: means after tide-iono-s1ab correction)')
        l1 = 'tide'
        l2 = 'final'
    else:
        l1 = 'tide'
        l2 = 'iono'
    export_esds2kml(framespd, esds, level1 = l1, level2 = l2, kmzfile = outkmzfile, overwrite = True, clean = False)
    print('done')

#%% main
if __name__ == "__main__":
    sys.exit(main())


