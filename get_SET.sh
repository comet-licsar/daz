
#this is to be run in bash, takes long (last full processing took 2 days)

# this is to get earth tides (bash):
in_esds=esds2021_take4.txt
in_frames=framespd_2021.csv
out_SET=earthtides_take4.csv

echo "frame, epoch, dEtide, dNtide, dUtide" > $out_SET
for aline in `cat $in_frames | tail -n+2 `; do
  frame=`echo $aline | cut -d ',' -f1`
  if [ `grep -c $frame $in_esds` -gt 0 ]; then
     echo $frame
     lon=`echo $aline | cut -d ',' -f3`
     lat=`echo $aline | cut -d ',' -f4`
     masterdate=`echo $aline | cut -d ',' -f2`
     masterdate=`echo ${masterdate:0:4}-${masterdate:4:2}-${masterdate:6:2}`
     centertime=`echo $aline | cut -d ',' -f9`
     masterdt=$masterdate"T"$centertime
     heading=`echo $aline | cut -d ',' -f5`
     mtide=`gmt earthtide -L$lon/$lat -T$masterdt 2>/dev/null | sed 's/\t/,/g'`
     NM=`echo $mtide | cut -d ',' -f2`
     EM=`echo $mtide | cut -d ',' -f3`
     VM=`echo $mtide | cut -d ',' -f4`
     
     for eline in `grep ^$frame $in_esds`; do
      #epochdate=`echo $eline | cut -d ',' -f2`
      epochdate=`echo $eline | cut -d ',' -f3`
      epochdatee=`echo ${epochdate:0:4}-${epochdate:4:2}-${epochdate:6:2}`
      epochdt=$epochdatee"T"$centertime
      etide=`gmt earthtide -L$lon/$lat -T$epochdt 2>/dev/null | sed 's/\t/,/g'`
      NE=`echo $etide | cut -d ',' -f2`
      EE=`echo $etide | cut -d ',' -f3`
      VE=`echo $etide | cut -d ',' -f4`
      U=`python3 -c "print(("$VE")-("$VM"))"`
      E=`python3 -c "print(("$EE")-("$EM"))"`
      N=`python3 -c "print(("$NE")-("$NM"))"`
      echo $frame","$epochdate","$E","$N","$U >> $out_SET
     done
  fi
done
