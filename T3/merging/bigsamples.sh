#threshold=4718592 # 4.5 GB
threshold=15728640 # 15 GB

if [ "$1" == "" ]
then
  samples=`cat $SUBMIT_WORKDIR/list_all.cfg | sed 's/_[[:digit:]]\+ .*//g'|sort -u`
else
  samples=$1
fi
bigsamples=""
smallsamples=""
for sample in $samples 
do
  size=`du -c $SUBMIT_OUTDIR/${sample}_*.root|grep total|sed 's/\t\+total//g'`
  if [ "$size" -gt "$threshold" ]
  then
    bigsamples="$bigsamples $sample"
  else
    smallsamples="$smallsamples $sample"
  fi
done
echo "BIG SAMPLES: $bigsamples"
echo "SMALL SAMPLES: $smallsamples"
