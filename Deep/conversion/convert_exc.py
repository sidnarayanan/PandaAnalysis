#!/usr/bin/env python

from sys import argv, exit
import numpy as np 
from os import getenv, system
from PandaCore.Utils.logging import logger
import PandaAnalysis.Deep.job_deepgen_utilities as deep_utils
from glob import glob 

fcfg = open(argv[1])
name = argv[2]
singletons = None
outdir = getenv('SUBMIT_NPY')
datadir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/deep/'
submit_name = getenv('SUBMIT_NAME')
me = argv[0].split('/')[-1]
argv = []
system('mkdir -p tmp/')


data = {}
for fpath in fcfg.readlines():
    print fpath.strip()
    try:
        d = np.load(fpath.strip())
    except [IOError, AttributeError] as e:
        logger.error(me, str(e))
        continue
    for k,v in d.iteritems():
        if k == 'singleton_branches':
            data[k] = v 
            continue
        if v.shape[0]:
            if k not in data:
                data[k] = []
            data[k].append(v)

if not len(data):
    logger.info(me, 'This was an empty config!')
    exit(0)

for k,v in data.iteritems():
    if k == 'singleton_branches':
        continue
    data[k] = np.concatenate(v)


def dump():
    global singletons 

    outpath = 'tmp/' + name + '_%s.npy'

    # particles
    d = data['particles']
    np.save(outpath%'particles', d)


dump()

ret = None
for ftmp in glob('tmp/*npy'):
    cmd = 'cp -v %s %s/train'%(ftmp,outdir)
    # cmd = 'gfal-copy -f file://$PWD/%s srm://t3serv006.mit.edu:8443/srm/v2/server?SFN=%s/%s'%(ftmp,outdir,ftmp.replace('tmp/',''))
    logger.info(me, cmd)
    ret = max(ret, system(cmd))

system('rm -rf tmp')
logger.debug(me, 'exit code %i'%ret)
exit(ret)
