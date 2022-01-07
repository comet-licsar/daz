# daz

Scripts to assess Sentinel-1 azimuth shift values to correct for non-displacements and decompose displacements from azimuth to N-S/E-W, originally developed to assess LiCSAR produced data.  

## daz_01_prepare_input.sh

Script to generate esds.txt and frames.txt files in LiCSAR environment.
This is an internal procedure partly covered by prepare_daz_licsar.sh that is to be exchanged for extraction of offsets directly from LiCSInfo database.

## daz_02_extract_SET.sh

Script to use GMT EarthTide (implementation of 'solid') to extract SET-related azimuth offsets.  
Output is esds.csv and frames.csv.

## daz_03_extract_iono.py

Approach to use IRI2016 model to extract and store STEC and Hiono and apply them to correct for daz_iono.

## daz_04_extract_PMM.py

Use of UNAVCO-hosted ITRF2014 PMM model to extract plate motion -related daz.

## daz_05_calculate_slopes.py

Script to use given column(s) (e.g. daz_mm, daz_mm_notide_noiono) to calculate azimuth velocity.

## daz_06_decompose.py

Script to decompose azimuth shifts into E,N motion components, optionally extract ITRF2014 PMM mean model values for given grid cell.

## daz_export2kmz.py

Script to export time series plots of daz values to KMZ file.  

For figures of decomposed data etc. see the daz_plotting library.


for binder see:
https://mybinder.org/v2/gl/comet_licsar%2Fdaz/HEAD
