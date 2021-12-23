#!/usr/bin/env python3

from daz_lib import *

# see get_SET.sh

# now (finally) load esds and tides to python, merge etc.
framescsv = 'framespd_2021.csv'
esdscsv = 'esds2021_take4.txt'
tidescsv = 'earthtides_take4.csv'

earthtides = pd.read_csv(tidescsv)
esds = pd.read_csv(esdscsv)
framespd = pd.read_csv(framescsv)

esds = merge_tides(esds, framespd, earthtides)


# to back up before continuing:
#framespd.to_csv('backup_0623_framespd.csv')
#esds.to_csv('backup_0623_esds.csv')
