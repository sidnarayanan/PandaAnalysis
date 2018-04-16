#!/usr/bin/env python

import ROOT as root
import numpy as np

from sys import path, stderr

import PandaCore.Tools.root_interface as r
import json 
from os import getenv, environ
from PandaCore.Tools.Misc import *
environ['KERAS_BACKEND'] = 'tensorflow'
from keras.models import Model, load_model
import PandaAnalysis.T3.job_utilities as utils
from glob import glob

cmssw_base = getenv('CMSSW_BASE')
sname = 'Deep.job_deepgen_utilities'

SAVE = False
STORE = False
NORM = False


def tree_to_arrays(infilepath, treename='inputs'):
    f = root.TFile(infilepath)
    t = f.Get(treename)
    data = {}
    singleton_branches = [x.GetName() for x in t.GetListOfBranches()
                                       if (x.GetName() != 'kinematics')]
    singleton_branches = sorted(list(set(singleton_branches)))
    singletons = r.read_tree(t, branches=singleton_branches)
    for k in singleton_branches:
        data[k] = singletons[k]

    data['singleton_branches'] = np.array(singleton_branches)

    arr = r.read_tree(t, branches=['kinematics'])
    data['particles'] = np.array([x[0].tolist() for x in arr])

    return data 


def normalize_arrays(data, cat, infilepath=None):
    if NORM and data['pt'].shape[0]:
        payload = json.load(open(cmssw_base + '/src/PandaAnalysis/data/deep/normalization_%s.json'%cat))
        mu = []; sigma = []
        for x in payload:
            mu_ = x['mean']
            sigma_ = x['stdev']
            if sigma_ == 0:
                sigma_ = 1
            mu.append(mu_)
            sigma.append(sigma_)
        mu = np.array(mu, dtype=np.float32); sigma = np.array(sigma, dtype=np.float32)

        data[cat] -= mu 
        data[cat] /= sigma

    if SAVE and (infilepath is not None):
        np.savez(infilepath.replace('.root','.npz'), **data)
    

def run_model(infilepattern, outfilepath):
    N = len(glob(infilepattern.replace('%i','*')))
    predictions = []
    for i in xrange(N):
        infilepath = infilepattern % i
        data = tree_to_arrays(infilepath)
        normalize_arrays(data, 'particles', infilepath)
        utils.print_time('preprocessing '+infilepath)
        if not STORE:
            utils.cleanup(infilepath)
