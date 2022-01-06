#!/usr/bin/env python3
"""
v1.0 2022-01-03 Milan Lazecky, Leeds Uni

This script will extract ionosphere shifts using IRI2016 model.

===============
Input & output files
===============
Inputs :
 - frames.csv - contains data with heading:
frame,master,center_lon,center_lat,heading,azimuth_resolution,avg_incidence_angle,centre_range_m,centre_time,dfDC
 - esds.csv - contains data with heading:
,frame,orbits_precision,daz_tide_mm,epochdate,daz_mm,years_since_beginning,daz_mm_notide

Outputs :
 - esds_with_iono.csv - added iono columns
 - frames_with_iono.csv - added iono columns

=====
Usage
=====
daz_03_extract_iono.py [--indaz esds.csv] [--infra frames.csv] [--outfra frames_with_iono.csv] [--outdaz esds_with_iono.csv]


"""
#%% Change log
'''
v1.0 2022-01-03 Milan Lazecky, Uni of Leeds
 - Original implementation - based on codes from 2021-06-24
'''
from daz_lib import *
from daz_iono import *

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
    indazfile = 'esds.csv'
    inframesfile = 'frames.csv'
    outdazfile = 'esds_with_iono.csv'
    outframesfile = 'frames_with_iono.csv'
    
    #%% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help", "indaz =", "infra =", "outdaz =", "outfra ="])
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
            elif o == "--outdaz":
                outdazfile = a
            elif o == "--outfra":
                outframesfile = a
        
        if os.path.exists(outdazfile):
            raise Usage('output esds csv file already exists. Cancelling')
        if os.path.exists(outframesfile):
            raise Usage('output frames csv file already exists. Cancelling')
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
    print('extra data cleaning step - perhaps should add to another step (first?)')
    esds, framespd = df_preprepare_esds(esds, framespd, firstdate = '', countlimit = 25)
    print('performing the iono calculation')
    esds, framespd = extract_iono_full(esds, framespd)
    print('saving files')
    esds.to_csv(outdazfile)
    framespd.to_csv(outframesfile)
    print('done')

#%% main
if __name__ == "__main__":
    sys.exit(main())

