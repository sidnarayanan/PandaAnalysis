#!/bin/bash

: ${SUBMIT_WORKDIR:?"Need to set SUBMIT_WORKDIR"}
: ${SUBMIT_OUTDIR:?"Need to set SUBMIT_OUTDIR"}
: ${SUBMIT_LOGDIR:?"Need to set SUBMIT_LOGDIR"}

myName=$(echo $0 | sed -e "s?.*/??g")

PInfo -n "$myName" "Cleaning up staging areas..."
WD=$SUBMIT_WORKDIR
rm -rf $WD/*
rm -rf $SUBMIT_LOGDIR/*

doTar=0
filesetSize=20
while getopts ":tn:" opt; do
  case $opt in
    t)
      doTar=1
      ;;
    n)
      filesetSize=$OPTARG
      ;;
    :)
      PError -n "$myName"  "Option -n must specify number of files"
      exit 1
      ;;
  esac 
done

PInfo -n "$myName" "Acquiring configuration..."
wget -O ${WD}/list.cfg $PANDA_CFG
${CMSSW_BASE}/src/PandaAnalysis/T3/bin/configBuilder.py --infile ${WD}/list.cfg --outfile ${WD}/local.cfg --nfiles $filesetSize
cp -v ${WD}/list.cfg ${WD}/list_all.cfg 
cp -v ${WD}/local.cfg ${WD}/local_all.cfg 

cd $CMSSW_BASE/
if [[ $doTar == 1 ]]; then
  PInfo -n "$myName" "Tarring up CMSSW..."
  tar --exclude-vcs -chzf cmssw.tgz src python biglib bin lib objs test external # h = --dereference symlinks
  mv -v cmssw.tgz ${WD}
fi

PInfo -n "$myName" "Creating executable..."
cd ${CMSSW_BASE}/src/PandaAnalysis/T3/inputs/
cp -v ${SUBMIT_TMPL} ${WD}/skim.py
chmod 775 ${WD}/skim.py

PInfo -n "$myName" "Finalizing work area..."
voms-proxy-init -voms cms --valid 168:00
cp -v /tmp/x509up_u$UID $WD/x509up

cp -v ${CMSSW_BASE}/src/PandaAnalysis/T3/inputs/exec.sh ${WD}

PInfo -n "$myName" "Taking a snapshot of work area..."
cp -rvT ${WD} $SUBMIT_OUTDIR/workdir/

PInfo -n "$myName" "Done!"

# input files for submission: cmssw.tgz, skim.py, x509up, local.cfg. exec.sh is the executable
