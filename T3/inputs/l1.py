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
import PandaAnalysis.T3.job_utilities as utils
from PandaAnalysis.Flat.analysis import *

Load('PandaAnalysisFlat')
data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'
utils.MAXCOPY = 1

def fn(input_name, isData, full_path):
    
    logger.info(sname+'.fn','Starting to process '+input_name)
    # now we instantiate and configure the analyzer
    a = analysis("l1", verbose=False)
    a.inpath = input_name
    a.outpath = utils.input_to_output(input_name)
    a.datapath = data_dir
    a.isData = isData
    utils.set_year(a, 2016)

    skimmer = root.pa.L1Analyzer(a)

    return utils.run_Analyzer(skimmer, isData, a.outpath)


if __name__ == "__main__":
    sample_list = cb.read_sample_config('local.cfg',as_dict=False)
    to_run = None 
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
    
    wd = utils.isolate()
    utils.main(to_run, processed, fn)

    utils.hadd(processed.keys())
    utils.print_time('hadd')

    ret = utils.stageout(outdir,outfilename)
    utils.cleanup('*.root')
    utils.un_isolate(wd)
    utils.print_time('stageout and cleanup')
    if not ret:
        utils.report_done(lockdir,outfilename,processed)
        utils.cleanup('*.lock')
        utils.print_time('create lock')
    else:
        exit(-1*ret)

    exit(0)
