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
#cats = ['inputs', 'inputs_dr', 'inputs_etaphi', 'targets', 'misc']
cats = ['inputs', 'targets', 'misc']
for cat in cats:
    branches[cat] = []
    f = open(basedir+cat+'.cfg')
    for l in f.readlines():
        x = l.strip()
        if x.startswith('#'):
            continue
        if x.endswith('_[]'):
            branches[cat] += [x.replace('[]','%i'%i) for i in xrange(5)]
        else:
            branches[cat].append(x)

cut = 'jotFlav == 5 && fabs(jotEta)<2.5 && jotGenPt>10 && jotGenEta>-10'
all_branches = [cut]
for _,v in branches.iteritems():
    all_branches += v

xarr = read_files(filenames = [infile],
                  branches = all_branches)

# flatten them
flattened = {}
for b in all_branches:
    flattened[b] = np.concatenate(xarr[b])

mask = flattened[cut].astype(bool)

# now merge some
data = {}
for cat in cats:
    if cat.startswith('input'):
        data[cat] = np.vstack([flattened[b][mask] for b in branches[cat]]).T
        data['shape'+cat.replace('inputs','')] = data[cat].shape 
for b in branches['targets'] + branches['misc']:
    data[b] = flattened[b][mask]


fout = infile.split('/')[-1].replace('.root','.npz')
np.savez(fout, **data) 
system('mv -v '+fout+' '+outpath)
