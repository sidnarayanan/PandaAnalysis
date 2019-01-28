#!/usr/bin/env python

from PandaCore.Tools.script import * 
from PandaCore.Tools.root_interface import Selector 
import PandaCore.Tools.Functions
import numpy as np
import json 
import os

args = parse('--out', '--name', '--json')

try:
    os.makedirs(args.out)
except OSError:
    pass

with open(args.json) as jsonfile:
    payload = json.load(jsonfile)
    weight = payload['weight']
    basedir = payload['base']
    features = payload['features']
    more_features = payload['more_features']
    cut = payload['cut']
    for i,s in enumerate(payload['samples']):
        if s['name'] == args.name:
            samples = s['samples']
            y = i
            break
    else:
        logger.error(sys.argv[0], 'Could not identify process '+args.name)
        sys.exit(1)

s = Selector()
chain = root.TChain('events')
for sample in samples:
    chain.AddFile(basedir + '/' + sample + '.root')

logger.info(sys.argv[0], 'Reading files for process '+args.name)
branches = [weight, 'hbbm_dreg']
branches += more_features
branches += [(f, -10, 6) for f in features]
s.read_tree(chain, branches=branches, cut=cut)

N, M = s['hbbjtidx'].shape

hbbidx = np.zeros_like(s[features[0]])
for k in xrange(M):
    hbbidx[np.arange(N), s['hbbjtidx'][:,k]] = 1

# temporary, will eventually be jotPt-dimensional
breg = np.ones_like(hbbidx)
for k in xrange(s['jotBReg'].shape[1]):
    breg[np.arange(N), s['hbbjtidx'][:,k]] = s['jotBReg'][:,k]

selected = []
for f in features:
    arr = np.empty((N,M))
    for k in xrange(M):
        arr[:,k] = s[f][np.arange(N),s['hbbjtidx'][:,k]]
    selected.append(arr)

#X = np.array([s[f] for f in features] + [hbbidx, breg]).transpose((1,2,0))

m = np.einsum('i,ij->ij', s['hbbm_dreg'], np.ones_like(s['jotBReg']))

X = np.array(selected + [s['jotBReg'], m]).transpose((1,2,0))

m = s['hbbm_dreg'].reshape((m.shape[0],1))
X = np.array([m]).transpose((1,2,0))

X = X.reshape((X.shape[0], X.shape[1]*X.shape[2]))
print X.shape
W = s[weight]
#W *= 1000 / W.sum()
#Y = y * np.ones(shape=W.shape)
Y = s['hbbm_dreg']

def save(arr, label):
    fout = args.out+'/'+args.name+'_'+label+'.npy'
    np.save(fout, arr)
    logger.info(sys.argv[0], 'Saved to '+fout)

save(X, 'x')
save(Y, 'y')
save(W, 'w')
