#!/usr/bin/env python 

from sys import argv, stdout, exit
from os import getenv, system 
from glob import glob 
import numpy as np 
from PandaAnalysis.Deep.NH1 import NH1
import json 
from pprint import pprint

limit = int(argv[1]) if len(argv) > 1 else None
datadir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/deep/'
outdir = getenv('SUBMIT_OUTDIR')

argv = []
import ROOT as root 
f_pt_scaled = root.TFile.Open(datadir + 'flatten_gen_scaled.root')
h_pt_scaled = f_pt_scaled.Get('h_QCD')

files = glob(outdir + '/QCD*.npz')
np.random.shuffle(files)
if limit:
    files = files[:limit]

NBINS = 100 

def reweight_s(x_pt):
    return h_pt_scaled.GetBinContent(h_pt_scaled.FindBin(x_pt))
reweight_s = np.vectorize(reweight_s)

branches = { }
bins = { }
cats = ['particles']
for cat in cats:
    print 'category',cat
    maxs = None; mins = None
    print '  acquiring ranges...'
    for i,f in enumerate(files[:20]):
        print '  %i/20\r'%i,; stdout.flush()
        arr = np.load(f)[cat]
        N = arr.shape[-1]
        tmp_maxs = np.array([np.percentile(arr[:,:,x].flatten(),98) for x in xrange(N)])
        tmp_mins = np.array([np.percentile(arr[:,:,x].flatten(),2) for x in xrange(N)])
        if maxs is None:
            maxs = tmp_maxs
            mins = tmp_mins
        else:
            maxs = np.maximum(tmp_maxs, maxs)
            mins = np.minimum(tmp_mins, mins)

    hists = []
    for i,(lo,hi) in enumerate(zip(mins,maxs)):
        if lo == hi: #???
            print 'idx %i has nearly constant output'%i
            hi = lo + 1
        width = float(hi-lo)/NBINS
        hists.append( NH1(np.arange(lo,hi+width,width)) )

    print '  filling histograms...'
    for i,f in enumerate(files):
        print '  %i/%i\r'%(i,len(files)),; stdout.flush()
        data = np.load(f)
        pt = data['pt']
        weight = reweight_s(pt)
        arr = data[cat]
        for ix in xrange(N):
            xarr = arr[:,:,ix]
            warr = np.array([weight for _ in xrange(xarr.shape[1])])
            xarr = xarr.flatten()
            warr = warr.flatten()
            hists[ix].fill_array(xarr, warr)
    
    dumps = []
    for ix in xrange(N):
        hist = hists[ix]
        if ix > 3:
            dumps.append({
                'idx'    : ix,
                'mean'   : hist.mean(),
                'stdev'  : hist.stdev(sheppard=True),
                'median' : hist.median(),
                'width'  : 0.5 * (hist.quantile(0.68) - hist.quantile(0.32))
            })
            if hist.stdev(sheppard=True) == 0:
                pprint(dumps[-1])
                print hist
        else: # kinematics, don't mess with them
            dumps.append({
                'idx'    : ix,
                'mean'   : 0,
                'stdev'  : 1,
                'median' : 0,
                'width'  : 1
            })

    json.dump(dumps, open(datadir + '/normalization_%s.json'%cat,'w'), indent=2)

