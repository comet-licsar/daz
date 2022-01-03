#!/bin/bash
# 2021-06-24
# M Lazecky 2020, 2021
# to generate all ESD values within frames..

#parameter is number of process, second from and third to track no

# to generate the ESDs - just export from esd table in LiCSInfo.
# this script is on how to create esds txt file in the form of:
#
# frame,esd_master,epoch,daz_total_wrt_orbits,daz_cc_wrt_orbits,orbits_precision,version
# 011A_04903_171717,20191004,20190407,0.0092,0.01161,P,m
# ...

# it is ok to have it just:
# frame,esd_master,epoch,daz_total_wrt_orbits,daz_cc_wrt_orbits

if [ -z $3 ]; then 
  echo "Please provide parameters: ID relorb_from relorb_to"; 
  echo "e.g. prepare_daz_licsar.sh 2021_all 1 175"
  exit;
fi

outesd=$LiCSAR_procdir/esds_$1.txt
outfr=$LiCSAR_procdir/esds_frames_$1.txt


#extra file, not used further for daz analysis:
outdif=$LiCSAR_procdir/esds_diff_in_frames_$1.txt




echo '' > $outdif
echo "frame,esd_master,epoch,daz_total_wrt_orbits,daz_cc_wrt_orbits,orbits_precision,version" > $outesd
echo "frame,master,center_lon,center_lat" > $outfr

# to run for a subset of tracks (used for parallel extraction)
for tr in `seq $2 $3`; do
#for tr in `seq 1 175`; do
for f in `ls $tr`; do
echo $f
m=`ls $tr/$f/SLC | head -n1`

if [ `ls $tr/$f/log/co*qu*log 2>/dev/null | wc -l` -gt 0 ]; then
#use epochs that are still existing / not qual-removed
rm $LiCSAR_procdir/$tr/$f/epochs.list 2>/dev/null
if [ ! `ls $LiCSAR_procdir/$tr/$f/IFG | wc -l` -eq `ls $LiCSAR_public/$tr/$f/interferograms | wc -l` ]; then 
 echo "number of geo vs ifg differs"
 echo $f >> $outdif
fi
for x in `ls $LiCSAR_procdir/$tr/$f/IFG | sed 's/_/ /'`; do echo $x >> $LiCSAR_procdir/$tr/$f/epochst.list; done
for x in `ls $LiCSAR_public/$tr/$f/interferograms | sed 's/_/ /'`; do echo $x >> $LiCSAR_procdir/$tr/$f/epochst.list; done
sort -u $LiCSAR_procdir/$tr/$f/epochst.list > $LiCSAR_procdir/$tr/$f/epochs.list
rm $LiCSAR_procdir/$tr/$f/epochst.list

hgtfile=$LiCSAR_public/$tr/$f/metadata/$f'.geo.hgt.tif'
ll=`gdalinfo $hgtfile | grep ^Center`
lon=`echo $ll | cut -d "," -f1 | cut -d '(' -f2`
lat=`echo $ll | cut -d ")" -f1 | cut -d ',' -f2`
echo $f","$m","$lon","$lat >> $outfr

for c in `ls $tr/$f/log/co*qu*log`; do
cc=''
esdtz=''
s=`basename $c | cut -d '_' -f4 | cut -d '.' -f1`
if [ `grep -c $s $LiCSAR_procdir/$tr/$f/epochs.list` -gt 0 ]; then
modc=`stat -c %y $c | gawk {'print $1'} | sed 's/-//g'`
if [ `datediff $s $modc` -lt 21 ]; then orb=R; else orb=P; fi
SM=$m
if [ `grep -c "estimation will be performed through " $c` -gt 0 ]; then
 SM=`grep "estimation will be performed through " $c | rev | gawk {'print $1'} | rev`
fi
esdtz=`grep "Total azimuth offset" $c | cut -d ':' -f2 | gawk {'print $1'} | cut -c -8`
match_it=0
#one version of GAMMA gives '_'.... 
keyword='matching_iteration_'
match_it=`grep -c $keyword $c`
let no=$match_it/2   
if [ $match_it -eq 0 ]; then
 #another version without '_'
 keyword='matching iteration '
 match_it=`grep -c "$keyword" $c`
 no=$match_it
fi
if [ $match_it -gt 0 ]; then
 #hmm, so the new version shows TOTAL VALUE ONLY FOR ESD towards the CC_intens-ied rslc...
 #let no=$match_it/2
 cmdadd=''
 for ma in `seq 1 $no`; do
  if [ ! -z `echo $keyword | cut -c 9` ]; then
   cc=`grep "$keyword"$ma $c  | gawk {'print $2'}`
  else
   cc=`grep -A1 "$keyword"$ma $c | tail -n1 | gawk {'print $3'}`
  fi
  cmdadd=$cmdadd"+("$cc")"
 done
 cc=`python -c "print(round("$cmdadd",6))" 2>/dev/null`
 esdtz=`python -c "print(round("$esdtz$cmdadd",6))" 2>/dev/null`
 ver=m
else
 cc=`grep "intensity_matching" $c | head -n1 | cut -d ':' -f2 | gawk {'print $1'}`
 ver=i
fi


if [ ! -z $esdtz ]; then
 echo $f","$SM","$s","$esdtz","$cc","$orb","$ver >> $outesd
fi
fi

done #coreg_quality log
fi #number of coreg qual logs
done #frame
done #track
