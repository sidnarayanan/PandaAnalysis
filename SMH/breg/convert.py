#!/usr/bin/env python

from sys import argv
from os import getenv, system
me = argv[0]
infile = argv[1]
outpath = getenv('SUBMIT_NPY')
basedir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/SMH/breg/'
argv = []

from PandaCore.Tools.root_interface import read_files, rename_dtypes
import PandaCore.Tools.Functions
import numpy as np

branches = {}
for cat in ['inputs', 'targets', 'misc']:
    branches[cat] = []
    f = open(basedir+cat+'.cfg')
    for l in f.readlines():
        x = l.strip()
        if x.endswith('_[]'):
            branches[cat] += [x.replace('[]','%i'%i) for i in xrange(5)]
        else:
            branches[cat].append(x)


all_branches = []
for _,v in branches.iteritems():
    all_branches += v

xarr = read_files(filenames = [infile],
                  branches = all_branches,
                  cut = 'jotFlav == 5 && fabs(jotEta)<2.5')

# flatten them
flattened = {}
for b in all_branches:
    flattened[b] = np.concatenate(xarr[b])

# now merge some
data = {}
data['inputs'] = np.vstack([flattened[b] for b in branches['inputs']]).T
for b in branches['targets'] + branches['misc']:
    data[b] = flattened[b]
data['shape'] = data['inputs'].shape 


np.savez(infile.split('/')[-1].replace('.root','.npz'), **data) 
system('mv -v *npz '+outpath)
