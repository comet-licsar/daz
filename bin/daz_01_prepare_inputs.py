#!/usr/bin/env python3

from daz_lib import *
from daz_lib_licsar import *


#######################################
# step 1 - export ESDs to a txt file:
# see prepare_daz_licsar.sh
# or just get it from licsinfo's esd table..
# (see LiCSquery - get_daz etc.)

dazfile = 'esds.txt'

#######################################
# step 2 - add additional information to the esds.txt and frames.txt files
# to be run at JASMIN
# this is to prepare the framespd and esds

fname = 'esds2021_frames.txt'
outcsv = 'framespd_2021.csv'
if not os.path.exists(outcsv):
    generate_framespd(fname, outcsv)

