#!/usr/bin/env python

from re import sub
from sys import argv,exit
from os import system,getenv,path
from time import clock,time
import json
from glob import glob

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
import PandaAnalysis.Deep.job_deep_utilities as deep_utils
from PandaAnalysis.Flat.analysis import deep

deep_utils.STORE = True
deep_utils.SAVE = True
deep_utils.INFER = False
deep_utils.NORM = False # temporary, need to recalculate normalizations

Load('PandaAnalyzer')
data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'

class BDTAdder(object):
    def __init__(self):
        Load('TMVABranchAdder')
        self.ba = root.TMVABranchAdder()
        self.ba.defaultValue = -1.2
        self.ba.presel = 'fj1ECFN_2_4_20>0'
        for v in tagcfg.variables:
            self.ba.AddVariable(v[0],v[2])
        for v in tagcfg.formulae:
            self.ba.AddFormula(v[0],v[2])
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
    skimmer = root.PandaAnalyzer()

    processType = utils.classify_sample(full_path, isData)
    if processType == root.kSignal:
        processType = root.kTop
    analysis = deep() 
    analysis.processType=processType 
    analysis.dump()
    skimmer.SetAnalysis(analysis)
    skimmer.isData=isData
    skimmer.AddPresel(root.FatJet450Sel())

    outpath = utils.run_PandaAnalyzer(skimmer, isData, input_name)
    if not outpath:
        return False 
    for f in glob('*_pf_*.root'):
        add_bdt(f, 'inputs')
    deep_utils.run_model(outpath.replace('.root','_pf_%i.root'), outpath)
    return True


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
    
    utils.main(to_run, processed, fn)

    utils.hadd(processed.keys())
    if deep_utils.STORE and False:
        utils.hadd([x.replace('output_','') for x in glob('*pf*.root')], 'arrays.root')
        utils.cleanup('*pf*.root')
    utils.print_time('hadd')

    add_bdt()
    utils.print_time('bdt')

    # utils.record_inputs('output.root',processed)
    # utils.print_time('record inputs')

    ret = utils.stageout(outdir,outfilename)
    if deep_utils.STORE and False:
        utils.stageout(outdir,outfilename.replace('.root','_arrays.root'),'arrays.root')
    utils.cleanup('*.root')
    if deep_utils.SAVE:
        data = {}
        for f in glob('*npz'):
            f_data = deep_utils.np.load(f)
            for k,v in f_data.iteritems():
                if k not in data:
                    data[k] = []
                if v.shape[0] > 0:
                    data[k].append(v)
        if len(data['pt']) > 0:
            merged_data = {k : deep_utils.np.concatenate(v) for k,v in data.iteritems()}
            deep_utils.np.savez('merged_arrays.npz', **merged_data)
            utils.print_time('merging npz')
            ret = max(ret, utils.stageout(outdir, outfilename.replace('.root', '.npz'), 'merged_arrays.npz'))
            #    utils.stageout(outdir,outfilename.replace('.root','.npz'),'arrays.npz')
        utils.cleanup('*.npz')
    utils.print_time('stageout and cleanup')
    if not ret:
        utils.report_done(lockdir,outfilename,processed)
        utils.cleanup('*.lock')
        utils.print_time('create lock')
    else:
        exit(-1*ret)

    exit(0)
