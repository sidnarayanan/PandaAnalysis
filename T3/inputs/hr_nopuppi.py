#!/usr/bin/env python

from os import getenv

from PandaCore.Utils.root import root 
from PandaCore.Tools.Misc import *
from PandaCore.Utils.load import *
import PandaAnalysis.T3.job_utilities as utils
from PandaAnalysis.Flat.analysis import * 

Load('PandaAnalyzer')
data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'

class BDTAdder(object):
    def __init__(self):
        Load('TMVABranchAdder')
        self.ba = root.TMVABranchAdder()
        self.ba.defaultValue = -1.2
        self.ba.presel = 'clf_ECFN_2_4_20>0'
        for v in utils.tagcfg.variables:
            self.ba.AddVariable(v[0],v[2].replace('fj','clf_').replace('[0]',''))
        for v in utils.tagcfg.formulae:
            self.ba.AddFormula(v[0],v[2].replace('fj','clf_').replace('[0]',''))
        for s in utils.tagcfg.spectators:
            self.ba.AddSpectator(s[0])
        self.ba.BookMVA('clf_top_ecf_bdt',data_dir+'/trainings/top_ecfbdt_v8_BDT.weights.xml')
    def __call__(self, fname='output.root', tname='events'):
        # now run the BDT
        self.ba.treename = tname
        self.ba.RunFile(fname)

def fn(input_name, isData, full_path):
    a = analysis('substructure',
                 recalcECF=False,
                 reclusterFJ=True,
                 puppiJets=False,
                 ak=False) 
    a.inpath = input_name
    a.outpath = utils.input_to_output(input_name)
    a.datapath = data_dir
    a.isData = isData
    utils.set_year(a, 2016)
    a.processType = utils.classify_sample(full_path, isData)	

    skimmer = root.pa.HRAnalyzer(a)

    return utils.run_HRAnalyzer(skimmer, isData, a.outpath)


if __name__ == "__main__":
    utils.wrapper(fn, post_fn=BDTAdder())

