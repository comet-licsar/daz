#!/usr/bin/env python3
"""
v1.0 2022-01-03 Milan Lazecky, Leeds Uni

This script will calculate solid Earth tides (currently using 'solid' implemented in GMT earthtides) and merge with esds.txt.
If the SET file exists, it will just merge it.

===============
Input & output files
===============
Inputs :
 - frames.csv - contains data with heading:
frame,master,center_lon,center_lat,heading,azimuth_resolution,avg_incidence_angle,centre_range_m,centre_time,dfDC
 - esds.txt - contains data with heading:
frame,esd_master,epoch,daz_total_wrt_orbits,daz_cc_wrt_orbits,orbits_precision,version
 - [ optional, if exists ] tides.csv - contains data with heading:
frame, epoch, dEtide, dNtide, dUtide

Outputs :
 - esds.csv - contains data with heading:
,frame,orbits_precision,daz_tide_mm,epochdate,daz_mm,years_since_beginning,daz_mm_notide
 - tides.csv - contains data with heading:
frame, epoch, dEtide, dNtide, dUtide

=====
Usage
=====
daz_02_extract_SET.py [--indaz esds.txt] [--infra frames.csv] [--tidescsv tides.csv] [--outdaz esds.csv]

 --tidescsv - input or output (if does not exist) file containing SET.

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
    indazfile = 'esds.txt'
    inframesfile = 'frames.csv'
    outdazfile = 'esds.csv'
    tidescsv = 'earthtides.csv'
    
    #%% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help", "indaz=", "infra=", "outdaz=", "tidescsv="])
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
            elif o == "--tidescsv":
                tidescsv = a
        
        if os.path.exists(outdazfile):
            raise Usage('output esds csv file already exists. Cancelling')
        if not os.path.exists(inframesfile):
            raise Usage('input frames txt file does not exist. Cancelling')
        if not os.path.exists(indazfile):
            raise Usage('input esds txt file does not exist. Cancelling')
            
    except Usage as err:
        print("\nERROR:",)
        print("  "+str(err.msg))
        print("\nFor help, use -h or --help.\n")
        return 2
    
    # processing itself:
    if not os.path.exists(tidescsv):
        print('SET file {0} does not exist. Generating it.'.format(tidescsv))
        print('(warning - this may take really long. it can take days..)')
        cmd = 'get_SET.sh {0} {1} {2}'.format(indazfile, inframesfile, tidescsv)
        os.system(cmd)
    else:
        print('SET file already exists. Will use it for merging')
    
    if not os.path.exists(tidescsv):
        print('ERROR - the SET file was not generated, exiting')
        return 2
    
    # now (finally) load esds and tides to python, merge etc.
    earthtides = pd.read_csv(tidescsv)
    esds = pd.read_csv(indazfile)
    framespd = pd.read_csv(inframesfile)
    print('converting SET data to azimuth direction and merging with ESD values')
    esds = merge_tides(esds, framespd, earthtides)
    print('exporting final merge to '+outdazfile)
    esds.to_csv(outdazfile)
    print('done')

#%% main
if __name__ == "__main__":
    sys.exit(main())
