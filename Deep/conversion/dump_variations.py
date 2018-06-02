#!/usr/bin/env python

import numpy as np
from os import environ

basedir = environ['SUBMIT_NPY'] + '/train/'

d = []
for fname in ['QCD_0', 'Top_lo_0']:
    data = np.load(basedir + '/' + fname + '_singletons.npy')
    d.append(data)

data = np.concatenate(d)

n_var = data.shape[1]
props = []
for i in xrange(n_var):
    x = data[:,i]
    props.append((np.mean(x), np.std(x)))


print '['
for p in props:
    print '    '+str(p)+','
print ']'
