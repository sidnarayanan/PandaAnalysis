#!/bin/bash

WD=$PWD

env
python -c "import sys; import socket; sys.stderr.write('hostname = '+socket.gethostname()+'\n');"
uname -a 1>&2 
lsb_release -a 1>&2 
hostname 1>&2

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

ls
mv local*cfg local.cfg

#export X509_USER_PROXY=${PWD}/x509up
export HOME=$WD

RELEASE=$CMSSW_VERSION
scram p CMSSW $RELEASE
tar xzf cmssw.tgz -C $RELEASE

cd $RELEASE
eval `scram runtime -sh`
cd -

echo -n "file length "
wc -l local.cfg

python skim.py $@

ls
rm -rf $RELEASE skim.py x509up cmssw.tgz local.cfg *root *npz lcg-cp* testfile
