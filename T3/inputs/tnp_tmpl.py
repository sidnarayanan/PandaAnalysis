#!/usr/bin/env python

from re import sub
from sys import argv,exit
from os import system,getenv,path
from time import clock,time
import json

which = int(argv[1])
submit_id = int(argv[2])
sname = argv[0]
argv=[]

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Utils.load import *
import PandaCore.Tools.job_config as cb
import PandaAnalysis.Tagging.cfg_v8 as tagcfg
import PandaAnalysis.T3.job_utilities as utils
from PandaAnalysis.Flat.analysis import *

Load('PandaAnalyzer')
data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'

class BDTAdder(object):
    def __init__(self):
        Load('TMVABranchAdder')
        self.ba = root.TMVABranchAdder()
        self.ba.defaultValue = -1.2
        self.ba.presel = 'fjECFN_2_4_20>0'
        for v in tagcfg.variables:
            self.ba.AddVariable(v[0],v[2].replace('fj1','fj'))
        for v in tagcfg.formulae:
            self.ba.AddFormula(v[0],v[2].replace('fj1','fj'))
        for s in tagcfg.spectators:
            self.ba.AddSpectator(s[0])
        self.ba.BookMVA('top_ecf_bdt',data_dir+'/trainings/top_ecfbdt_v8_BDT.weights.xml')
    def __call__(self, fname='output.root', tname='events'):
        # now run the BDT
        self.ba.treename = tname
        self.ba.RunFile(fname)

add_bdt = BDTAdder() #backwards compatability


def fn(input_name, isData, full_path):
    
    logger.info(sname+'.fn','Starting to process '+input_name)
    # now we instantiate and configure the analyzer
    a = monotop(True)
    a.recalcECF = True 
    a.varyJESTotal = True
    a.mcTriggers = True
    a.inpath = input_name
    a.outpath = utils.input_to_output(input_name)
    a.datapath = data_dir
    a.isData = isData
    utils.set_year(a, 2016)
    a.processType = utils.classify_sample(full_path, isData)	

    skimmer = root.pa.PandaAnalyzer(a)
    skimmer.AddPresel(root.pa.FatJetSel())

    return utils.run_PandaAnalyzer(skimmer, isData, a.outpath)


if __name__ == "__main__":
    sample_list = cb.read_sample_config('local.cfg',as_dict=False)
    to_run = None #sample_list[which]
    for s in sample_list:
        if which==s.get_id():
            to_run = s
            break
    if not to_run:
        logger.error(sname,'Could not find a job for PROCID=%i'%(which))
        exit(3)
    
    outdir = getenv('SUBMIT_OUTDIR')
    lockdir = getenv('SUBMIT_LOCKDIR')  
    outfilename = to_run.name+'_%i.root'%(submit_id)
    processed = {}

    utils.report_start(outdir,outfilename,to_run.files)
    
    utils.main(to_run, processed, fn)

    utils.hadd(processed.keys())
    utils.print_time('hadd')

    add_bdt()
    utils.print_time('bdt')

    ret = utils.stageout(outdir,outfilename)
    utils.cleanup('*.root')
    utils.print_time('stageout and cleanup')
    if not ret:
        utils.report_done(lockdir,outfilename,processed)
        utils.cleanup('*.lock')
        utils.print_time('create lock')
    else:
        exit(-1*ret)

    exit(0)
