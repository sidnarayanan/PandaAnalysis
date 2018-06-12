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
import PandaCore.Tools.job_management as cb
import PandaAnalysis.Tagging.cfg_v8 as tagcfg

Load('PandaAnalyzer')
data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'

stopwatch = clock() 
def print_time(label):
    global stopwatch
    now_ = clock()
    logger.debug(sname+'.print_time:'+str(time()),
           '%.3f s elapsed performing "%s"'%((now_-stopwatch)/1000.,label))
    stopwatch = now_

def copy_local(long_name):
    replacements = {
                r'\${EOS}':'root://eoscms.cern.ch//store/user/snarayan',
                r'\${EOS2}':'root://eoscms.cern.ch//store/group/phys_exotica',
                r'\${CERNBOX}':'root://eosuser//eos/user/s/snarayan',
                r'\${CERNBOXB}':'root://eosuser//eos/user/b/bmaier',
            }
    full_path = long_name
    for k,v in replacements.iteritems():
        full_path = sub(k,v,full_path)
    logger.info(sname,full_path)

    panda_id = long_name.split('/')[-1].split('_')[-1].replace('.root','')
    input_name = 'input_%s.root'%panda_id

    # if the file is cached locally, why not use it?
    local_path = full_path.replace('root://xrootd.cmsaf.mit.edu/','/mnt/hadoop/cms')
    if path.isfile(local_path):
        # apparently SmartCached files can be corrupted...
        ftest = root.TFile(local_path)
        if ftest and not(ftest.IsZombie()):
            full_path = local_path

    '''
    # xrdcp if remote, copy if local - DEPRECATED
    if 'root://' in full_path:
        system('xrdcp %s %s'%(full_path,input_name))
    else:
        system('cp %s %s'%(full_path,input_name))
    '''

    # rely on pxrdcp for local and remote copies
    # default behavior: drop PF candidates
    cmd = "pxrdcp %s %s '!pfCandidates'"%(full_path,input_name)
    logger.info(sname+'.copy_local',cmd)

    system(cmd)
            
    if path.isfile(input_name):
        logger.info(sname+'.copy_local','Successfully copied to %s'%(input_name))
        return input_name
    else:
        logger.error(sname+'.copy_local','Failed to copy %s'%input_name)
        return None


def fn(input_name,isData,full_path):
    
    logger.info(sname+'.fn','Starting to process '+input_name)
    # now we instantiate and configure the analyzer
    skimmer = root.PandaAnalyzer()
    skimmer.isData=isData
    skimmer.SetFlag('firstGen',True)
    skimmer.AddPresel(root.MonotopSel())
    processType=root.PandaAnalyzer.kNone
    if not isData:
        if any([x in full_path for x in ['Vector_','Scalar_']]):
            processType=root.PandaAnalyzer.kSignal
        elif any([x in full_path for x in ['ST_','ZprimeToTT']]):
            processType=root.PandaAnalyzer.kTop
        elif 'EWKZ2Jets' in full_path:
            processType=root.PandaAnalyzer.kZEWK
        elif 'EWKW' in full_path:
            processType=root.PandaAnalyzer.kWEWK
        elif 'ZJets' in full_path or 'DY' in full_path:
            processType=root.PandaAnalyzer.kZ
        elif 'WJets' in full_path:
            processType=root.PandaAnalyzer.kW
        elif 'GJets' in full_path:
            processType=root.PandaAnalyzer.kA
        elif 'TTJets' in full_path or 'TT_' in full_path:
            processType=root.PandaAnalyzer.kTT
    skimmer.processType=processType 

    # read the inputs
    try:
        fin = root.TFile.Open(input_name)
        tree = fin.FindObjectAny("events")
        weight_table = fin.FindObjectAny('weights')
        hweights = fin.FindObjectAny("hSumW")
    except:
        logger.error(sname+'.fn','Could not read %s'%input_name)
        return False # file open error => xrootd?
    if not tree:
        logger.error(sname+'.fn','Could not recover tree in %s'%input_name)
        return False
    if not hweights:
        logger.error(sname+'.fn','Could not recover hweights in %s'%input_name)
        return False
    if not weight_table:
        weight_table = None

    output_name = input_name.replace('input','output')
    skimmer.SetDataDir(data_dir)
    if isData:
        with open(data_dir+'/certs/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt') as jsonFile:
            payload = json.load(jsonFile)
            for run_str,lumis in payload.iteritems():
                run = int(run_str)
                for l in lumis:
                    skimmer.AddGoodLumiRange(run,l[0],l[1])
    rinit = skimmer.Init(tree,hweights,weight_table)
    if rinit:
        logger.error(sname+'.fn','Failed to initialize %s!'%(input_name))
        return False 
    skimmer.SetOutputFile(output_name)

    # run and save output
    skimmer.Run()
    skimmer.Terminate()

    ret = path.isfile(output_name)
    if ret:
        logger.info(sname+'.fn','Successfully created %s'%(output_name))
        return True
    else:
        logger.error(sname+'.fn','Failed in creating %s!'%(output_name))
        return False


def cleanup(fname):
    ret = system('rm -f %s'%(fname))
    if ret:
        logger.error(sname+'.cleanup','Removal of %s exited with code %i'%(fname,ret))
    else:
        logger.info(sname+'.cleanup','Removed '+fname)
    return ret


def hadd(good_inputs):
    good_outputs = ' '.join([x.replace('input','output') for x in good_inputs])
    cmd = 'hadd -f output.root ' + good_outputs
    ret = system(cmd)    
    if not ret:
        logger.info(sname+'.hadd','Merging exited with code %i'%ret)
    else:
        logger.error(sname+'.hadd','Merging exited with code %i'%ret)


def add_bdt():
    # now run the BDT
    Load('TMVABranchAdder')
    ba = root.TMVABranchAdder()
    ba.treename = 'events'
    ba.defaultValue = -1.2
    ba.presel = 'fj1ECFN_2_4_20>0'
    for v in tagcfg.variables:
        ba.AddVariable(v[0],v[2])
    for v in tagcfg.formulae:
        ba.AddFormula(v[0],v[2])
    for s in tagcfg.spectators:
        ba.AddSpectator(s[0])
    ba.BookMVA('top_ecf_bdt',data_dir+'/trainings/top_ecfbdt_v8_BDT.weights.xml')
    ba.RunFile('output.root')


def drop_branches(to_drop=None, to_keep=None):
    # remove any irrelevant branches from the final tree.
    # this MUST be the last step before stageout or you 
    # run the risk of breaking something
    
    if not to_drop and not to_keep:
        return 0

    if to_drop and to_keep:
        logger.error(sname+'.drop_branches','Can only provide to_drop OR to_keep')
        return 0

    f = root.TFile('output.root','UPDATE')
    t = f.FindObjectAny('events')
    n_entries = t.GetEntriesFast() # to check the file wasn't corrupted
    if to_drop:
        if type(to_drop)==str:
            t.SetBranchStatus(to_drop,False)
        else:
            for b in to_drop:
                t.SetBranchStatus(b,False)
    elif to_keep:
        t.SetBranchStatus('*',False)
        if type(to_keep)==str:
            t.SetBranchStatus(to_keep,True)
        else:
            for b in to_keep:
                t.SetBranchStatus(b,True)
    t_clone = t.CloneTree()
    f.WriteTObject(t_clone,'events','overwrite')
    f.Close()

    # check that the write went okay
    f = root.TFile('output.root')
    if f.IsZombie():
        logger.error(sname+'.drop_branches','Corrupted file trying to drop '+to_drop)
        return 1 
    t_clone = f.FindObjectAny('events')
    if (n_entries==t_clone.GetEntriesFast()):
        return 0
    else:
        logger.error(sname+'.drop_branches','Corrupted tree trying to drop '+to_drop)
        return 2


def stageout(outdir,outfilename):
    mvargs = 'mv $PWD/output.root %s/%s'%(outdir,outfilename)
    logger.info(sname,mvargs)
    ret = system(mvargs)
    system('rm *.root')
    if not ret:
        logger.info(sname+'.stageout','Move exited with code %i'%ret)
    else:
        logger.error(sname+'.stageout','Move exited with code %i'%ret)
        return ret
    if not path.isfile('%s/%s'%(outdir,outfilename)):
        logger.error(sname+'.stageout','Output file is missing!')
        ret = 1
    return ret


def report_done(outdir,outfilename,processed):
    flock = open(outdir+'/locks/'+outfilename.replace('.root','.lock'),'w')
    for k,v in processed.iteritems():
        flock.write(v+'\n')
    flock.close()


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

    outdir = 'XXXX' # will be replaced when building the job
    outfilename = to_run.name+'_%i.root'%(submit_id)
    processed = {}
    
    print_time('loading')
    for f in to_run.files:
        input_name = copy_local(f)
        print_time('copy %s'%input_name)
        if input_name:
            success = fn(input_name,(to_run.dtype!='MC'),f)
            print_time('analyze %s'%input_name)
            if success:
                processed[input_name] = f
            cleanup(input_name)
            print_time('remove %s'%input_name)
    
    if len(processed)==0:
        exit(1)

    hadd(list(processed))
    print_time('hadd')

    add_bdt()
    print_time('add BDT')
    if drop_branches(to_drop='fj1ECFN*'):
        exit(2)
    print_time('drop branches')
    ret = stageout(outdir,outfilename)
    print_time('stageout')
    if not ret:
        report_done(outdir,outfilename,processed)
        print_time('create lock')
    else:
        exit(-1*ret)

    exit(0)
