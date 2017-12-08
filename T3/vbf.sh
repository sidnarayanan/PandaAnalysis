#!/bin/bash

export PATH=${PATH}:${CMSSW_BASE}/src/PandaCore/bin/

export PANDA="${CMSSW_BASE}/src/PandaAnalysis"
export PANDA_CFG="http://t3serv001.mit.edu/~snarayan/histcatalog/20171129_lo.cfg" 
#export PANDA_CFG="http://t3serv001.mit.edu/~snarayan/histcatalog/20170912_vbfH.cfg" 
#export PANDA_CFG="http://t3serv001.mit.edu/~snarayan/histcatalog/20170830_vbf.cfg" 
export PANDA_FLATDIR="/data/t3home000/zdemirag/forSid/panda/v_004_12"
mkdir -p $PANDA_FLATDIR

export PANDA_VBFSCAN=${PANDA_FLATDIR}/scan/
mkdir -p $PANDA_VBFSCAN

export SUBMIT_TMPL="skim_vbf_tmpl.py"
export SUBMIT_NAME="v_004_12"
export SUBMIT_WORKDIR="/data/t3home000/zdemirag/forSid/condor/"${SUBMIT_NAME}"/work/"
export SUBMIT_LOGDIR="/data/t3home000/zdemirag/forSid/condor/"${SUBMIT_NAME}"/logs/"
export SUBMIT_OUTDIR="/data/t3home000/zdemirag/forSid/store/"${SUBMIT_NAME}"/batch/"
mkdir -p $SUBMIT_WORKDIR $SUBMIT_OUTDIR/locks/ $SUBMIT_LOGDIR

export SUBMIT_LOCALACCESS=1
