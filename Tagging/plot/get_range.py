#!/usr/bin/env python

from PandaCore.Tools.script import * 
from PandaCore.Tools.root_interface import Selector 
from PandaAnalysis.Monotop.CombinedBVetoSelection import cuts 
import numpy as np 
import json 

betas = [0.5, 1., 2., 4.]
orders = [1, 2, 3]
Ns = [2, 3, 4]

f = root.TFile.Open( getenv('PANDA_FLATDIR') + '/WJets.root' )
t = f.Get('events')

class ECF(object):
    def __init__(self, b, o, N):
        self.b = b
        self.o = o
        self.N = N
    def __str__(self):
        return 'fjECFN_%i_%i_%.2i[0]'%(self.o, self.N, int(10*self.b))
    def label(self):
        return 'e(%i,%i,%.1f)'%(self.o, self.N, self.b)

class ECFR(object):
    def __init__(self, n, d):
        self.n = n
        self.d = d 
        self.x = (n.b*n.o)/(d.b*d.o)
    def eval(self, selector):
        if np.min(selector[str(self.d)]) == 0:
            return np.zeros(1)
        return (selector[str(self.n)]) / np.power(selector[str(self.d)], self.x)
    def formula(self):
        return '%s/pow(%s, %.2f)'%(str(self.n), str(self.d), self.x)
    def label(self):
        return '%s/(%s)^{%.2f}'%(self.n.label(), self.d.label(), self.x)
    def __repr__(self):
        return 'lambda a: a["%s"]/np.power(a["%s"],%f)'%(str(self.n), str(self.d), self.x)


ecfs = [ECF(b, o, N) for b in betas
                     for o in orders
                     for N in Ns]
s = Selector()
print 'Reading...'
s.read_tree(t, branches=[str(e) for e in ecfs],
            cut=tAND(cuts['singlemuonw'], 'fjPt[0]>250 && fjECFN_2_4_20[0]>0'))
print '...Read'

ecf_ranges = [{'formula':str(e), 'label':e.label(),
               'lo':np.min(s[str(e)]), 'hi':np.max(s[str(e)])} for e in ecfs]

#for e in ecfs:
#    print e.label(), np.min(s[str(e)]), np.max(s[str(e)])

ratios = {}

for n in ecfs:
    for d in ecfs:
        if d.N > n.N:
            continue
        r = ECFR(n, d)
        vals = r.eval(s)
        std = vals.std()
        hi = vals.max()
        hi = min(hi, vals.mean()+(5*std))
        ratios[r] = {
                'formula' : repr(r),
                'label' : r.label(),
                'lo' : np.min(vals),
                'hi' : hi,
                'median' : np.percentile(vals, 50),
            }

with open('ranges.json','w') as fp:
    json.dump({'ratios':ratios.values(), 'ecfs':ecf_ranges}, 
               fp, sort_keys=True, indent=4, separators=(',',':'))
