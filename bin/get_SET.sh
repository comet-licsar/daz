#!/bin/bash
# 2023-08-31: Muhammet Nergizci: improved vartype to keep using bc (faster than through calling python)
# 2021-06-24
# M Lazecky 2020, 2021
# this is to get earth tides (bash):
# NOTE: takes long..
#
# output csv has epochdate column as YYYYMMDD

if [ -z $3 ]; then
 echo "Usage: in_esds in_frames out_SET"
 echo "e.g.: esds.txt frames.csv tides.csv"
 exit
fi


#in_esds=esds2021_take4.txt
#in_frames=framespd_2021.csv
#out_SET=earthtides_take4.csv

in_esds=$1
in_frames=$2
out_SET=$3

if [ ! -f $in_esds ] || [ ! -f $in_frames ] || [ -f $out_SET ]; then
  echo "ERROR. Please check if following files exist:"
  echo "in_esds: "$in_esds
  echo "in_frames: "$in_frames
  echo "note this file must NOT exist:"
  echo "out_SET: "$out_SET
  exit
fi

j='none'
# checking first the epochtime column:
getetime=0
i=0; for x in `head -n1 $in_esds | sed 's/\,/ /g'`; do let i=$i+1; if [ $x == 'epochtime' ]; then j=$i; fi; done
if [ $j == 'none' ]; then
  # if not, revert to original 'inaccurate' SET estimation
  i=0; for x in `head -n1 $in_esds | sed 's/\,/ /g'`; do let i=$i+1; if [ `echo $x | cut -c -5` == 'epoch' ]; then j=$i; fi; done
else
  getetime=1
fi
if [ $j == 'none' ]; then
 echo "ERROR, input file does not contain epoch or epochdate column"
 exit
fi

echo "frame,epoch,dEtide,dNtide,dUtide" > $out_SET
for aline in `cat $in_frames | tail -n+2 `; do
  frame=`echo $aline | cut -d ',' -f1`
  if [ `grep -c $frame $in_esds` -gt 0 ]; then
     echo $frame
     lon=`echo $aline | cut -d ',' -f3`
     lat=`echo $aline | cut -d ',' -f4`
     masterdate=`echo $aline | cut -d ',' -f2`
     masterdate=`echo ${masterdate:0:4}-${masterdate:4:2}-${masterdate:6:2}`
     #if [ $getetime == 0 ]; then
     centertime=`echo $aline | cut -d ',' -f9`
     masterdt=$masterdate"T"$centertime
     #else
     #  masterdt=`grep $masterdate $in_esds | cut -d ',' -f $j | sed 's/ /T/'`
     #  echo "debug - assuming this as masterdt"
     #  echo $masterdt
     #fi
     heading=`echo $aline | cut -d ',' -f5`
     mtide=`gmt earthtide -L$lon/$lat -T$masterdt 2>/dev/null | sed 's/\t/,/g'`
     NM=`echo $mtide | cut -d ',' -f2`
     EM=`echo $mtide | cut -d ',' -f3`
     VM=`echo $mtide | cut -d ',' -f4`
     ##Convert scientific notation to decimal notation (e.g. e^3 to 0.001)
     NM=$(printf "%.12f" "$NM")
     EM=$(printf "%.12f" "$EM")
     VM=$(printf "%.12f" "$VM")


     for eline in `grep ^$frame $in_esds | sed 's/ /T/'`; do
      #epochdate=`echo $eline | cut -d ',' -f2`
      if [ $getetime == 0 ]; then
        epochdate=`echo $eline | cut -d ',' -f$j`
        if [ `echo $epochdate | grep -c '-'` == 0 ]; then
          epochdatee=`echo ${epochdate:0:4}-${epochdate:4:2}-${epochdate:6:2}`
        else
          epochdatee=$epochdate
          epochdate=`echo $epochdate | sed 's/-//g'`
        fi
        epochdt=$epochdatee"T"$centertime
      else
        epochdate=`echo $eline | cut -d ',' -f$j | cut -d 'T' -f 1 | cut -d ' ' -f 1 | sed 's/-//g'` # in case there is space, although now it should not
        epochdt=`echo $eline | cut -d ',' -f$j | sed 's/ /T/'`
      fi
      etide=`gmt earthtide -L$lon/$lat -T$epochdt 2>/dev/null | sed 's/\t/,/g'`
      NE=`echo $etide | cut -d ',' -f2`
      EE=`echo $etide | cut -d ',' -f3`
      VE=`echo $etide | cut -d ',' -f4`
      
      ##Convert scientific notation to decimal notation
      NE=$(printf "%.12f" "$NE")
      EE=$(printf "%.12f" "$EE")
      VE=$(printf "%.12f" "$VE")
      
      U=`echo '('$VE')-('$VM')' | bc`
      E=`echo '('$EE')-('$EM')' | bc`
      N=`echo '('$NE')-('$NM')' | bc`
      #U=`python3 -c "print(("$VE")-("$VM"))"`
      #E=`python3 -c "print(("$EE")-("$EM"))"`
      #N=`python3 -c "print(("$NE")-("$NM"))"`
      echo $frame","$epochdate","$E","$N","$U >> $out_SET
     done
  fi
done
