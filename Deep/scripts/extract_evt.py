#!/usr/bin/env python

import numpy as np
from glob import glob
from sys import argv, exit
from PandaCore.Tools.Misc import PInfo, PWarning

flist = glob(argv[2])
target_evt = np.int64(argv[1])

lumis = []
pts = []
msds = []
pfs = []

for f in flist:
    arr = np.load(f)
    evt = arr['eventNumber']
    mask = (evt == target_evt)
    if np.sum(mask):
        idx = np.argmax(mask)
        pfs.append(arr['pf'][idx])
        pts.append(arr['pt'][idx])
        msds.append(arr['msd'][idx])
        lumis.append(arr['lumi'])
        PInfo(argv[0], 'Found %i in %s'%(target_evt, f))


if lumis:
    np.savez('sliced.npz', pf=np.array(pfs), msd=np.array(msds), pt=np.array(pts), lumi=np.array(lumis))
else:
    PError(argv[0], 'Could not find %i in %s'%(target_evt, argv[2]))
    exit(1)
