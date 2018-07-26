#!/usr/bin/env python

from sys import argv
from os import getenv, system
from re import sub
me = argv[0]
infile = argv[1]
outpath = getenv('SUBMIT_NPY')
basedir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/SMH/breg/'
argv = []

from PandaCore.Tools.root_interface import read_files, rename_dtypes
import PandaCore.Tools.Functions
import numpy as np

def add_idx(x, i):
    return sub(r'jot([A-z]*)',r'jot\1[hbbjtidx[%i]]'%i, x)

branches = {}
cats = ['inputs', 'targets', 'misc_higgs']
#cats = ['inputs', 'inputs_dr', 'inputs_etaphi', 'targets', 'misc']
for cat in cats:
    f = open(basedir+cat+'.cfg')
    for l in f.readlines():
        x = l.strip()
        if x.startswith('#'):
            continue
        if x.endswith('_[]'):
            branches.setdefault(cat+'_hbb0',[]).extend([add_idx(x.replace('[]','%i'%i), 0) for i in xrange(5)])
            branches.setdefault(cat+'_hbb1',[]).extend([add_idx(x.replace('[]','%i'%i), 1) for i in xrange(5)])
        elif x.endswith('[0]'):
            branches.setdefault(cat+'_hbb0',[]).append(x)
        elif x.endswith('[1]'):
            branches.setdefault(cat+'_hbb1',[]).append(x)
        else:
            branches.setdefault(cat+'_hbb0',[]).append(add_idx(x,0))
            branches.setdefault(cat+'_hbb1',[]).append(add_idx(x,1))

cut = 'hbbm>0'
all_branches = [cut]
for _,v in branches.iteritems():
    all_branches += v

xarr = read_files(filenames = [infile],
                  branches = all_branches,
                  cut = cut)

# flatten them
flattened = {}
for b in all_branches:
#    flattened[b] = np.concatenate(xarr[b])
    flattened[b] = xarr[b].astype(float)

# now merge some
data = {}
for cat in branches:
    if cat.startswith('input'):
        data[cat] = np.vstack([flattened[b] for b in branches[cat]]).T
        data['shape'+cat.replace('inputs','')] = data[cat].shape 
    else:
        for b in branches[cat]:
            data[b] = flattened[b]


fout = infile.split('/')[-1].replace('.root','.npz')
np.savez(fout, **data) 
system('mv -v '+fout+' '+outpath)
