#!/usr/bin/env python3
"""
v1.0 2022-01-03 Milan Lazecky, Leeds Uni

This script is to be run at JASMIN. It will grab output of prepare_daz_licsar.sh (esds.txt and frames.txt) or using LiCSquery - get_daz etc.,
extract additional data form the LiCSAR system and save as frames.csv (note esds.txt is not really needed now, as its csv will be prepared in next steps).

===============
Input & output files
===============
Inputs :
 - frames.txt - contains data with heading:
frame,master,center_lon,center_lat
 [not used now] esds.txt - contains data with heading:
frame,esd_master,epoch,daz_total_wrt_orbits,daz_cc_wrt_orbits,orbits_precision,version

Outputs :
 - frames.csv - contains data with heading:
frame,master,center_lon,center_lat,heading,azimuth_resolution,avg_incidence_angle,centre_range_m,centre_time,dfDC

=====
Usage
=====
daz_01_prepare_inputs.py [--infra frames.txt] [--outfra frames.csv]

 -...

"""
#%% Change log
'''
v1.0 2022-01-03 Milan Lazecky, Uni of Leeds
 - Original implementation - based on codes from 2021-06-24
'''

import getopt, os, sys

from daz_lib import *
try:
    from daz_lib_licsar import *
except:
    print('error importing licsar daz libraries. perhaps you are not at JASMIN?')
    exit()


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
    #indazfile = 'esds.txt'
    inframesfile = 'frames.txt'
    #outdazfile = 'esds.csv'
    outframesfile = 'frames.csv'
    
    #%% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help", "inesd =", "infra =", "outesd =", "outfra ="])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            #elif o == "--inesd":
            #    indazfile = a
            elif o == "--infra":
                inframesfile = a
            #elif o == "--outesd":
            #    outdazfile = a
            elif o == "--outfra":
                outframesfile = a

        if os.path.exists(outframesfile):
            raise Usage('output frames csv file already exists. Cancelling')
        if not os.path.exists(inframesfile):
            raise Usage('input frames txt file does not exist. Cancelling')

    except Usage as err:
        print("\nERROR:",)
        print("  "+str(err.msg))
        print("\nFor help, use -h or --help.\n")
        return 2
    
    # processing itself:
    generate_framespd(inframesfile, outframesfile)

#%% main
if __name__ == "__main__":
    sys.exit(main())
