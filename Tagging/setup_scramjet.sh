#!/bin/bash

export SCRAMJET=${CMSSW_BASE}/src/SCRAMJet
export TOP=${SCRAMJET}/TopTagging

#location of production 
histDir=scramjet
export SCRAMJETHIST=${EOS}/$histDir

#location of flat trees
#export SCRAMJETFLAT=${HOME}/home000/scramjet/v5/
#export SCRAMJETFLAT=${HOME}/scratch5/store/scramjet/v7/
export SCRAMJETFLAT=${HOME}/home000/store/scramjet/v7/
export SUBMIT_OUTDIR=${HOME}/hadoop/scramjet/v6/
