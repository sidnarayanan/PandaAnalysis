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
sname = 'Deep.job_deep_utilities'

SAVE = False
STORE = False
NORM = True
INFER = False

singleton_branches = ['msd','pt', 'rawpt', 'eta', 'phi',  'eventNumber',
                      'partonM', 'partonPt', 'partonEta', 'nPartons',
                      'nBPartons', 'nCPartons',
                      'rho','rawrho','rho2','rawrho2',
                      'tau32','tau32SD','tau21','tau21SD',
                      'maxcsv','mincsv','doubleb', # uncomment in a bit
                      'top_ecf_bdt']

def tree_to_arrays(infilepath, treename='inputs'):
    f = root.TFile(infilepath)
    t = f.Get(treename)
    data = {}
    singletons = r.read_tree(t, branches=singleton_branches)
    for k in singleton_branches:
        data[k] = singletons[k]

    arr = r.read_tree(t, branches=['kinematics'])
    data['pf'] = np.array([x[0].tolist() for x in arr])

    arr = r.read_tree(t, branches=['svs'])
    data['sv'] = np.array([x[0].tolist() for x in arr])

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
    

def infer(data):
    model = load_model(cmssw_base + '/src/PandaAnalysis/data/deep/regularized_conv_abel.h5')
    prediction = model.predict([data['pf'], data['msd'], data['pt']])
    prediction = np.array([prediction[:,2], prediction[:,3]]).T
    prediction = np.array(prediction, dtype = [('dnn_higgs',np.float32), ('dnn_top',np.float32)]) 
        # higgs, top quark
    return prediction


def arrays_to_tree(outfilepath, prediction, treename='dnn'):
    f = root.TFile(outfilepath, 'UPDATE')
    t = r.array_as_tree(prediction, treename=treename, fcontext=f)
    f.WriteTObject(t, treename, 'Overwrite')
    f.Close()


def run_model(infilepattern, outfilepath):
    N = len(glob(infilepattern.replace('%i','*')))
    predictions = []
    for i in xrange(N):
        infilepath = infilepattern % i
        data = tree_to_arrays(infilepath)
        normalize_arrays(data, 'pf', infilepath)
        normalize_arrays(data, 'sv', infilepath)
        utils.print_time('preprocessing')
        if INFER:
            pred = infer(data)
            if pred.shape[0]:
                predictions.append(pred)
            utils.print_time('inference')
        if not STORE:
            utils.cleanup(infilepath)
    if INFER:
        if predictions:
            pred = np.concatenate(predictions)
        else:
            pred = np.array([])
        arrays_to_tree(outfilepath, pred)
        utils.print_time('saving prediction')

