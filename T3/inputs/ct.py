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

Load('PandaAnalysisFlat')
data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'

def fn(input_name, isData, full_path):
    a = analysis('graph', ak=True, ak8=True, puppiJets=True)
    a.inpath = input_name
    a.outpath = utils.input_to_output(input_name)
    a.datapath = data_dir
    a.isData = isData
    utils.set_year(a, 2016)

    skimmer = root.pa.CTAnalyzer(a)

    return utils.run_Analyzer(skimmer, isData, a.outpath)

def post_fn():
    sleep(60)
    f = root.TFile.Open('output.root')
    t = f.Get('events')
    logger.debug('post_fn', 'Tree has %i entries'%(t.GetEntries()))

    logger.info('post_fn', 'Creating numpy arrays')
    arr = read_files(['output.root'], branches=None)
    adj = arr['adj']
    x = np.stack([arr['node'+f] for f in ['Pt','Eta','Phi','E','IsFinal','IsRoot']], axis=-1)
    jets = np.stack([arr['jet'+f] for f in ['Tau32','Tau21','MSD','Pt','Eta','Phi','M','PdgId']], axis=-1)
    logger.info('post_fn', 'Saving numpy arrays')
    np.savez_compressed('output.npz', adj=adj, x=x, jets=jets)
    do('mv -v output.npz output.root')

if __name__ == "__main__":
    utils.wrapper(fn, post_fn=post_fn) 

