#!/usr/bin/env python3
"""
v1.0 2022-01-03 Milan Lazecky, Leeds Uni

This script is to be run at JASMIN. It will grab output files (esds.txt and frames.txt) of either daz_lib_licsar.extract_all2txt() or (older) prepare_daz_licsar.sh script.
extract additional data form the LiCSAR system and save as frames.csv (note esds.txt is not really needed now, as its csv will be prepared in next steps).

===============
Input & output files
===============
Inputs :
 - frames.txt - contains data with heading:
frame,master,center_lon,center_lat
 - esds_orig.txt - contains data with heading:
frame,esd_master,epoch,daz_total_wrt_orbits,daz_cc_wrt_orbits,orbits_precision,version

Outputs :
 - frames.csv - contains data with heading:
frame,master,center_lon,center_lat,heading,azimuth_resolution,avg_incidence_angle,centre_range_m,centre_time,dfDC
 - esds.txt - same as before, but would be corrected for the shift between datasets pre/post 2020-07-29/30 if you used flag --orbdiff_fix

=====
Usage
=====
daz_01_prepare_inputs.py [--infra frames.txt] [--outfra frames.csv] [--indaz esds_orig.txt] [--outdaz esds.txt] [--orbdiff_fix]

 --orbdiff_fix - would apply fix due to change in orbits in 2020-07-29/30.  If working in LiCSAR environment, it will apply real difference, otherwise will apply 39 mm constant shift.
 (note the shift varies from this average by +-2std=25 mm and we observed also introduced bias in velocity e.g. 2 mm/year)
"""
#%% Change log
'''
v1.2 2023-08-30 ML
 - improved the orbit diff correction by using real difference between orbits
v1.1 2022-12-07 Milan Lazecky
 - include better correction of esds.txt related to the change in orbits. here we will use the 39 mm shift to align data before and after 2020-07-29 (or 30)
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
    indazfile = 'esds_orig.txt'
    inframesfile = 'frames.txt'
    outdazfile = 'esds.txt'
    outframesfile = 'frames.csv'
    orbdiff_fix = False
    
    #%% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help", "orbdiff_fix", "indaz=", "infra=", "outdaz=", "outfra="])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == "--orbdiff_fix":
                orbdiff_fix = True
                try:
                    from orbit_lib import *
                    using_orbits = True
                except:
                    print('WARNING: LiCSAR orbit library was not loaded. Fixing orbits using only constant value of 39 mm.')
                    using_orbits = False
            elif o == "--indaz":
                indazfile = a
            elif o == "--infra":
                inframesfile = a
            elif o == "--outdaz":
                outdazfile = a
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
    print('Generating output frames csv file - ETA couple of hours as we derive lot of info from frame txt files etc.')
    framespd = generate_framespd(inframesfile, outframesfile)
    
    # working with esds file
    esds, framespd = load_csvs(esdscsv = indazfile, framescsv = outframesfile)
    
    if orbdiff_fix:
        print('fixing the orb diff values in '+indazfile)
        esds = fix_pod_offset(esds, using_orbits=using_orbits)
    
    # some other details to prepare:
    try:
        esds = flag_s1b_esds(esds, framespd)
    except:
        print('unable to flag S1A/B for now, skipping')
    esds.to_csv(outdazfile, index=False)
    #else:
    #    # just reload it and save - at least will check for consistency
    #    esds=pd.read_csv(indazfile)
    #    esds.to_csv(outdazfile, index=False)


#%% main
if __name__ == "__main__":
    sys.exit(main())
