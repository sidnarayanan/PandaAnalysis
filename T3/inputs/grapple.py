#!/usr/bin/env python

from os import getenv

from PandaCore.Utils.root import root 
from PandaCore.Tools.Misc import *
from PandaCore.Utils.load import *
import PandaAnalysis.T3.job_utilities as utils
from PandaAnalysis.Flat.analysis import * 
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

    return utils.run_PandaAnalyzer(skimmer, isData, a.outpath)

def post_fn():
    do('hadd -f output.root ' + ' '.join(glob('*aux*root'))) 

if __name__ == "__main__":
    utils.wrapper(fn, post_fn=post_fn) 

