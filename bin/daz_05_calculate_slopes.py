#!/usr/bin/env python3
"""
v1.0 2022-01-03 Milan Lazecky, Leeds Uni

This script will estimate velocities from daz values.

===============
Input & output files
===============
Inputs :
 - frames.csv - should have ITRF and iono corr.
 - esds.csv - should have iono corr.

Outputs :
 - esds_final.csv
 - frames_final.csv

=====
Usage
=====
daz_05_calculate_slopes.py [--s1ab] [--indaz esds_with_iono.csv] [--infra frames_with_itrf.csv] [--outfra frames_final.csv] [--outdaz esds_final.csv]

Parameters:
    --s1ab ... also estimate (and store to outfra) the s1ab offset prior to velocity estimation. Now done only for the noiono+notide (final) daz
"""
#%% Change log
'''
v1.1 2024-04-06 ML
 - added S1AB offset estimation
v1.0 2022-01-03 Milan Lazecky, Uni of Leeds
 - Original implementation - based on codes from 2021-06-24
'''
from daz_lib import *
from daz_timeseries import *


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
    indazfile = 'esds_with_iono.csv'
    inframesfile = 'frames_with_itrf.csv'
    outdazfile = 'esds_final.csv'
    outframesfile = 'frames_final.csv'
    # for skipping roll assistance - as e.g. Turkey is not correct in ITR2014 PMM
    roll_assist = True
    s1ab = False
    
    #%% Read options
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "h", ["help", "s1ab", "indaz=", "infra=", "outdaz=", "outfra="])
        except getopt.error as msg:
            raise Usage(msg)
        for o, a in opts:
            if o == '-h' or o == '--help':
                print(__doc__)
                return 0
            elif o == "--s1ab":
                s1ab = True
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
    
    '''
    # step 6+ -- correct for daz_ARP=-39 mm: 29th July for S1A and 30th July for S1B
    #
    # update 2022-12-07 - now done in the step daz_01 already
    #################
    cols = ['daz_mm','daz_mm_notide', 'daz_mm_notide_noiono_grad']
    if ('s1AorB' in framespd.columns) and ('s1AorB' not in esds.columns):
        # ok, we can flag S1A/B and be more precise
        esds = flag_s1b_esds(esds, framespd)
        Bs = [esds['s1AorB']=='B']
        epB = Bs[Bs.epochdate => pd.Timestamp('2020-07-30')][cols]
        esds.update(epB.subtract(-39))
        As = [esds['s1AorB']=='A']
        epA = As[As.epochdate => pd.Timestamp('2020-07-29')][cols]
        esds.update(epA.subtract(-39))
    else:
        ep = esds[esds.epochdate => pd.Timestamp('2020-07-30')][cols]
        esds.update(ep.subtract(-39))
    '''
    # setting 'subset' - means, only data > 2016-03-01 as before it is too noisy
    #subset = True
    if subset:
        print('Subsetting dataset to include only data after 2016-03-01')
        esds = esds[esds['epochdate'] > pd.Timestamp('2016-03-01')]
    if s1ab:
        # estimate the offset first, then apply correction, and then use Huber as usual
        framespd = estimate_s1ab_allframes(esds, framespd, col = 'daz_mm_notide_noiono', rmsiter = 50)
        print('Applying S1AB corrections (only to daz_mm_notide_noiono and stored as daz_mm_final)')
        esds['daz_mm_final'] = esds['daz_mm_notide_noiono'].copy()
        esds, framespd = correct_s1ab(esds, framespd, cols=['daz_mm_final'])
    # 2021-10-12: the original way:
    for col in ['daz_mm', 'daz_mm_notide', 'daz_mm_notide_noiono_grad', 'daz_mm_notide_noiono_iri', 'daz_mm_notide_noiono','daz_mm_final']:
        if col in esds:
            print('estimating velocities of '+col)
            esds, framespd = df_calculate_slopes(esds, framespd, alpha = 1, eps = 1.35, bycol = col, subset = subset, roll_assist = roll_assist)
    #esds, framespd = df_calculate_slopes(esds, framespd, alpha = 1, eps = 1.35, bycol = 'daz_mm_notide_noiono_F2')
    # to back up before continuing:
    print('saving datasets')
    framespd.to_csv(outframesfile)
    esds.to_csv(outdazfile)
    print('done')

#%% main
if __name__ == "__main__":
    sys.exit(main())








