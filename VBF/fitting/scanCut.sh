#!/bin/bash

thiscut=$(echo $1 | sed "s?'??g")
scramdir=$2
outdir=$3
datacard=$4
label=$(echo $thiscut | sed "s?[\(\)_]??g" | sed "s?[><&]?_?g" )
thiscut=$(echo $thiscut | sed "s?&?\\\&?g")

echo "$thiscut"
echo "$scramdir" 
echo "$outdir"
echo "$label"

pwd
WD=$PWD

# cp -r $scramdir .
# ls
# cd CMSSW_7_4_7/src/MonoX
# scramv1 b ProjectRename
cd $scramdir/src/
# eval `scramv1 runtime -sh`
eval `/cvmfs/cms.cern.ch/common/scramv1 runtime -sh`
echo $CMSSW_BASE
which combine

cd $WD
cp -r $scramdir/src/MonoX .
cd MonoX
sed "s?XXXX?${thiscut}?g" configs/vbf_tmpl.py > configs/vbf_scan.py

python buildModel.py vbf_scan
python runModelNoPhoton.py

cp combined_model.root datacardsHiggs/ 
cd datacardsHiggs
#combine -M Asymptotic vbf_combined_fewbins.txt -n $label --run=blind --noFitAsimov 
#combine -M Asymptotic vbf_combined_binned.txt -n $label --run=blind --noFitAsimov 
#combine -M Asymptotic vbf_combined.txt -n $label --run=blind --noFitAsimov 
combine -M Asymptotic $datacard -n $label --run=blind --noFitAsimov 
ls higgs*root
mkdir -p $outdir
ret=$( cp higgs*root $outdir )

cd $WD
rm -rf MonoX 

exit $ret
