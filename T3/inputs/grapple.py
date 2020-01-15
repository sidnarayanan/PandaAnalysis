#!/usr/bin/env python

from os import getenv

from PandaCore.Utils.root import root 
from PandaCore.Tools.Misc import *
from PandaCore.Utils.load import *
import PandaAnalysis.T3.job_utilities as utils
from PandaAnalysis.Flat.analysis import * 
from PandaCore.Tools.root_interface import *
from PandaCore.Tools.script import do 
import numpy as np
from time import sleep
from glob import glob

Load('PandaAnalysisFlat')
data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'

def fn(input_name, isData, full_path):
    a = empty()
    a.inpath = input_name
    a.outpath = utils.input_to_output(input_name)
    a.datapath = data_dir
    a.isData = isData
    utils.set_year(a, 2017)
    a.processType = utils.classify_sample(full_path, isData)	

    skimmer = root.pa.GrappleAnalyzer(a)
    # skimmer.lastEvent = 10

    return utils.run_PandaAnalyzer(skimmer, isData, a.outpath)

def post_fn():
    do('hadd -f output.root ' + ' '.join(glob('*aux*root'))) 
    arr = read_files(['output.root'], branches=['kinematics', 'npv'],
                     treename='inputs')
    k = arr['kinematics']
    k = np.stack(k, axis=0)
    shape = k.shape
    k = np.stack(k.flatten(), axis=0)
    k = k.reshape(shape + (k.shape[-1],))
    # nmax = 200
    # evt = 3
    # print('k',k[1, :nmax, :])
    q = k[:,:,5]
    # print('q',q[evt, :nmax])
    y = k[:,:,4] == 0
    # print(y.astype(int).min())
    # print('y',y[evt, :nmax])
    x = k[:,:,:6]
    # print('pt',x[evt, :nmax, 0])
    # print('x',x[evt, :nmax, 4])
    x[:,:,4][(q == 0)] = -1 # if it's a neutral particle, mask the vertex ID
    # print('x',x[evt, :nmax, 4])
    p = k[:,:,-1]
    m = arr['genMET']
    np.savez('output.npz', x=x, y=y, q=q, p=p, m=m)
    do('mv -v output.npz output.root')

if __name__ == "__main__":
    utils.wrapper(fn, post_fn=post_fn) 

