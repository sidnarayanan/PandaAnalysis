#!/usr/bin/env python

import os
import json
import socket
import requests
from re import sub
from glob import glob 
from random import choice
from sys import exit, argv, stdout
from time import clock, time, sleep
from os import system, getenv, path, environ, getpid

from PandaCore.Tools.Misc import *
from PandaCore.Utils.load import Load
from PandaCore.Utils.root import root 
import PandaCore.Tools.job_config as cb
import PandaAnalysis.Tagging.cfg_v8 as tagcfg

_sname = 'T3.job_utilities'                                      # name of this module
_data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'    # data directory
_host = socket.gethostname()                                     # where we're running
_is_t3 = (_host.startswith('t3') and _host.endswith('mit.edu'))  # are we on the T3?
_remote_read = True                                              # should we read from hadoop or copy locally?
local_copy = bool(environ.get('SUBMIT_LOCALACCESS', True))       # should we always xrdcopy from T2?
_task_name = getenv('SUBMIT_NAME')                               # name of this task
_user = getenv('SUBMIT_USER')                                    # user running the task
_job_id = None                                                   # identifying string for this job
_year = 2016                                                     # what year's data is this analysis?
maxcopy = 3                                                      # maximum number of stagein attempts
_to_hdfs = bool(int(getenv('SUBMIT_HDFSCACHE', False)))          # should we cache on hdfs instead of local
_users = ['snarayan', 'bmaier', 'dhsu', 'ceballos', 'ballen']    # MIT T3 PandaAnalysis users 

stageout_protocol = None                                         # what stageout should we use?
if _is_t3:
    stageout_protocol = 'cp' 
else:
    # install certificates
    system('wget -nv http://t3serv001.mit.edu/~cmsprod/tgz/certificates.tar.gz')
    system('tar xvaf certificates.tar.gz')
    environ['X509_CERT_DIR'] = getenv('PWD') + '/certificates/'
    if system('which gfal-copy') == 0:
        stageout_protocol = 'gfal'
    elif True:
        logger.error(_sname, 
                    'Sorry, lcg is no longer supported - cannot trust T3 gridftp')
        raise RuntimeError 
    elif system('which lcg-cp') == 0:
        stageout_protocol = 'lcg'
    else:
        try:
            ret = system('wget -nv http://t3serv001.mit.edu/~snarayan/misc/lcg-cp.tar.gz')
            ret = max(ret, system('tar -xvf lcg-cp.tar.gz'))
            if ret:
                raise RuntimeError
            environ['PATH'] = '$PWD/lcg-cp:'+environ['PATH']
            environ['LD_LIBRARY_PATH'] = '$PWD/lcg-cp:'+environ['LD_LIBRARY_PATH']
            stageout_protocol = 'lcg'
        except Exception as e:
            logger.error(_sname, 
                         'Could not install lcg-cp in absence of other protocols!')
            raise e


# derived from t3serv006.mit.edu:/etc/bestman2/conf/bestman2.rc
_gsiftp_doors = [
        't3serv010.mit.edu', 
        ]



# global to keep track of how long things take
_stopwatch = time() 
def print_time(label):
    global _stopwatch
    now_ = time()
    logger.debug(_sname+'.print_time:', 
           '%.1f s elapsed performing "%s"'%((now_-_stopwatch), label))
    _stopwatch = now_

# set the data-taking period
def set_year(analysis, year):
    global _year
    analysis.year = year
    _year = year

# isolate the job
def isolate():
    pid = getpid()
    p = 'job_%i'%pid 
    try:
        os.mkdir(p)
    except OSError:
        pass
    os.chdir(p)
    return p

def un_isolate(p):
    os.chdir('..')
    cleanup(p)

# convert an input name to an output name
def input_to_output(name):
    if 'input' in name:
        return name.replace('input', 'output')
    else:
        return 'output_' + name.split('/')[-1]

def _validate_file(p):
    if not path.isfile(p):
        return False
    ftest = root.TFile(p)
    if bool(ftest) and not(ftest.IsZombie()):
        logger.info(_sname+'._validate_file', '%s is a good file'%p)
        return True 
    return False

def request_data(xrd_path, first_attempt):
    logger.info(_sname+'.request_data', xrd_path)

    panda_id = xrd_path.split('/')[-1].split('_')[-1].replace('.root', '')
    input_path = 'input_%s.root'%panda_id

    if 'scratch' in xrd_path:
        local_path = xrd_path.replace('root://t3serv006.mit.edu/', '/mnt/hadoop')
    else:
        local_path = xrd_path.replace('root://xrootd.cmsaf.mit.edu/', '/mnt/hadoop/cms')

    # if we're on the T3, only attempt to access the local data once, otherwise go for T2
    if first_attempt or not _is_t3:
        # first see if data is already present - if so, great!
        if _validate_file(local_path):
            logger.info(_sname+'.request_data', 'Using local file %s'%local_path)
            return local_path 

        if _is_t3:
            for user in _users:
                local_user_path = local_path.replace('paus', user)
                if _validate_file(local_user_path):
                    # tell server we're using this file
                    payload = {'path' : local_user_path, 
                               'bytes' : path.getsize(local_user_path)}
                    r = requests.post(cb.report_server+'/condor/requestdata', json=payload)
                    if r.status_code == 200:
                        logger.info(_sname+'.request_data', 'return=%s'%(str(r).strip()))
                    else:
                        logger.warning(_sname+'.request_data', 'return=%s'%(str(r).strip()))
                    logger.info(_sname+'.request_data', 'Using local user file %s'%local_path)
                    return local_user_path

    # ok now we have to copy the data in:
    xrdargs = ['xrdcopy', '-f', xrd_path]
    if not stdout.isatty():
        xrdargs.insert(1, '--nopbar')
    cache = _is_t3 and _to_hdfs
    if cache:
        input_path = local_path.replace('paus', _user) 
        parent = '/'.join(input_path.split('/')[:-1])
        if not path.isdir(parent):
            try:
                logger.info(_sname+'.request_data', 'creating parent at '+parent)
                os.makedirs(parent)
                os.chmod(parent, 0777)
            except OSError as e:
                logger.warning(_sname+'.request_data', str(e))
                pass 
    xrdargs.append(input_path)
    xrdargs = ' '.join(xrdargs)
    logger.info(_sname+'.request_data', xrdargs)
    ret = system(xrdargs)
    if ret:
        logger.error(_sname+'.request_data', 'Failed to xrdcopy %s'%input_path)
        return None 
    if _validate_file(input_path):
        if cache:
            payload = {'path' : input_path, 
                       'bytes' : path.getsize(input_path)}
            r = requests.post(cb.report_server+'/condor/requestdata', json=payload)
            os.chmod(input_path, 0777)
        logger.info(_sname+'.request_data', 'Successfully xrdcopied %s'%input_path)
        return input_path
    return None 


# wrapper around remove. be careful!
def cleanup(fname, _verbose=True):
    if path.isfile(fname):
        os.remove(fname)
    elif path.isdir(fname):
        os.rmdir(fname)
    else:
        for f in glob(fname):
            cleanup(f, False)
    if _verbose:
        logger.info(_sname+'.cleanup', 'Removed '+fname)
    return 0 # if it made it this far without OSError, it's good


# wrapper around hadd
def hadd(good_inputs, output='output.root'):
    good_outputs = ' '.join([input_to_output(x) for x in good_inputs])
    cmd = 'hadd -f ' + output + ' ' + good_outputs
    logger.info(_sname+'.hadd', cmd)
    ret = system(cmd)    
    if not ret:
        logger.info(_sname+'.hadd', 'Merging exited with code %i'%ret)
    else:
        logger.error(_sname+'.hadd', 'Merging exited with code %i'%ret)


# remove any irrelevant branches from the final tree.
# this MUST be the last step before stageout or you 
# run the risk of breaking something
def drop_branches(to_drop=None, to_keep=None):
    if not to_drop and not to_keep:
        return 0

    if to_drop and to_keep:
        logger.error(_sname+'.drop_branches', 'Can only provide to_drop OR to_keep')
        return 0

    f = root.TFile('output.root', 'UPDATE')
    t = f.FindObjectAny('events')
    n_entries = t.GetEntriesFast() # to check the file wasn't corrupted
    if to_drop:
        if type(to_drop)==str:
            t.SetBranchStatus(to_drop, False)
        else:
            for b in to_drop:
                t.SetBranchStatus(b, False)
    elif to_keep:
        t.SetBranchStatus('*', False)
        if type(to_keep)==str:
            t.SetBranchStatus(to_keep, True)
        else:
            for b in to_keep:
                t.SetBranchStatus(b, True)
    t_clone = t.CloneTree()
    f.WriteTObject(t_clone, 'events', 'overwrite')
    f.Close()

    # check that the write went okay
    f = root.TFile('output.root')
    if f.IsZombie():
        logger.error(_sname+'.drop_branches', 'Corrupted file trying to drop '+str(to_drop))
        return 1 
    t_clone = f.FindObjectAny('events')
    if (n_entries==t_clone.GetEntriesFast()):
        return 0
    else:
        logger.error(_sname+'.drop_branches', 'Corrupted tree trying to drop '+str(to_drop))
        return 2



# stageout a file (e.g. output or lock)
#  - if _is_t3, execute a simple cp
#  - else, use lcg-cp
# then, check if the file exists:
#  - if _is_t3, use os.path.isfile
#  - else, use lcg-ls
def stageout(outdir, outfilename, infilename='output.root', n_attempts=10):
    gsiftp_doors = _gsiftp_doors[:]
    if stageout_protocol is None:
        logger.error(_sname+'.stageout', 
               'Stageout protocol has not been satisfactorily determined! Cannot proceed.')
        return -2
    timeout = 300
    ret = -1
    for i_attempt in xrange(n_attempts):
#        door = choice(gsiftp_doors); gsiftp_doors.remove(door)
        door = gsiftp_doors[0]
        failed = False
        if stageout_protocol == 'cp':
            cpargs = ['cp', 
                      '-v', 
                      '$PWD/%s'%infilename, 
                      '%s/%s'%(outdir, outfilename)]
            lsargs = ['cp', cpargs[-1], '$PWD/testfile']
            rmargs = ['rm', cpargs[-1]]
        elif stageout_protocol == 'gfal':
            cpargs = ['gfal-copy', 
                      '-f', 
                      '--transfer-timeout %i'%timeout, 
                      '$PWD/%s'%infilename, 
                      'gsiftp://%s:2811//%s/%s'%(door, outdir, outfilename)]
            lsargs = ['gfal-stat', cpargs[-1]]
            rmargs = ['gfal-rm', cpargs[-1]]
        elif stageout_protocol == 'lcg':
            cpargs = ['lcg-cp', 
                      '-v -D srmv2 -b', 
                      'file://$PWD/%s'%infilename, 
                      'gsiftp://%s:2811//%s/%s'%(door, outdir, outfilename)]
            lsargs = ['lcg-ls', 
                      '-v -D srmv2 -b']

        cpargs = ' '.join(cpargs)
        lsargs = ' '.join(lsargs)
        rmargs = ' '.join(rmargs)

        logger.info(_sname+'.stageout', cpargs)
        ret = system(cpargs)
        if not ret:
            logger.info(_sname+'.stageout', 'Move exited with code %i'%ret)
            sleep(10) # give the filesystem a chance to respond
        else:
            logger.warning(_sname+'.stageout', 'Move exited with code %i'%ret)
            failed = True
        if not failed:
            logger.info(_sname+'.stageout', lsargs)
            ret = system(lsargs)
            if ret:
                logger.warning(_sname+'.stageout', 'Output file is missing!')
                failed = True
        if not failed:
            logger.info(_sname+'.stageout', 'Copy succeeded after %i attempts'%(i_attempt+1))
            cleanup('testfile')
            return ret
        else:
            system(rmargs)
            timeout = int(timeout * 1.5)
        cleanup('testfile')
    logger.error(_sname+'.stagoeut', 'Copy failed after %i attempts'%(n_attempts))
    return ret

# report home that the job has started
def report_start(outdir, outfilename, args):
    global _job_id
    if not cb.textlock:
        _job_id = '_'.join(outfilename.replace('.root', '').split('_')[-2:]) + '_' + str(int(time()))
        hashed = [cb.md5hash(x) for x in args]
        payload = {'starttime' : int(time()), 
                   'host' : _host, 
                   'task' : _task_name + '_' + _user, 
                   'job_id' : _job_id, 
                   'args' : hashed}
        for _ in xrange(5):
            try:
                r = requests.post(cb.report_server+'/condor/start', json=payload)
                if r.status_code == 200:
                    logger.info('T3.job_utilities.report_start', 'return=%s'%(str(r).strip()))
                    return 
                else:
                    logger.error('T3.job_utilities.report_start', 'return=%s'%(str(r).strip()))
                sleep(10)
            except requests.ConnectionError as e:
                logger.error('T3.job_utilities.report_start', str(e))
        logger.error('T3.job_utilities.report_start', 'Dying after 5 attempts!')
        raise requests.ConnectionError()


# write a lock file, based on what succeeded, 
# and then stage it out to a lock directory
def report_done(outdir, outfilename, processed):
    if cb.textlock:
        outfilename = outfilename.replace('.root', '.lock')
        flock = open(outfilename, 'w')
        for k, v in processed.iteritems():
            flock.write(v+'\n')
        flock.close()
        stageout(outdir, outfilename, outfilename)
        cleanup('*.lock')
    else:
        # this really has to work, so go forever
        while True:
            try:
                payload = {'timestamp' : int(time()), 
                           'task' : _task_name + '_' + _user, 
                           'job_id' : _job_id, 
                           'args' : [cb.md5hash(v) for _, v in processed.iteritems()]}
                logger.info('T3.job_utilities.report_done', 
                            'payload=\n%s'%repr(payload))
                r = requests.post(cb.report_server+'/condor/done', json=payload)
                if r.status_code == 200:
                    logger.info('T3.job_utilities.report_done', 'return=%s'%(str(r).strip()))
                    return 
                else:
                    logger.error('T3.job_utilities.report_done', 'return=%s'%(str(r).strip()))
                sleep(10)
            except requests.ConnectionError as e:
                logger.error('T3.job_utilities.report_done', str(e))


# make a record in the primary output of what
# inputs went into it
def record_inputs(outfilename, processed):
    fout = root.TFile.Open(outfilename, 'UPDATE')
    names = root.TNamed('record', 
                        ', '.join(processed.values()))
    fout.WriteTObject(names)
    fout.Close()


# classify a sample based on its name
def classify_sample(full_path, isData):
    _classification = [
                (root.pa.kSignal , ['Vector_', 'Scalar_']), 
                (root.pa.kTop    , ['ST_', 'ZprimeToTT']), 
                (root.pa.kZEWK   , 'EWKZ2Jets'), 
                (root.pa.kWEWK   , 'EWKW'), 
                (root.pa.kZ      , ['ZJets', 'DY', 'Z1Jets', 'Z2Jets']), 
                (root.pa.kW      , ['WJets', 'W1Jets', 'W2Jets']), 
                (root.pa.kA      , 'GJets'), 
                (root.pa.kTT     , ['TTJets', 'TT_', 'TTTo']), 
                (root.pa.kH      , 'HTo'), 
                (root.pa.kVV     , ['WW', 'WZ', 'ZZ', 'WpWp']), 
            ]
    if not isData:
        for e, pattern in _classification:
            if type(pattern) == str:
                if pattern in full_path:
                    return e 
            else:
                if any([x in full_path for x in pattern]):
                    return e
    return root.pa.kNoProcess

# add top-tagging BDT
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
        self.ba.BookMVA('top_ecf_bdt',_data_dir+'/trainings/top_ecfbdt_v8_BDT.weights.xml')
    def __call__(self, fname='output.root', tname='events'):
        # now run the BDT
        self.ba.treename = tname
        self.ba.RunFile(fname)



# read a CERT json and add it to the skimmer
_jsons = {
        2016 : '/certs/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt', 
        2017 : '/certs/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt', 
        2018 : '/certs/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt', 
        }
def add_json(skimmer):
    json_path = _jsons.get(_year, None)
    if not json_path:
        logger.error("T3.job_utilities.add_json", "Unknown key = "+str(_year))
    else:
        logger.info("T3.job_utilities.add_json", json_path)
    json_path = _data_dir + json_path
    with open(json_path) as jsonFile:
        payload = json.load(jsonFile)
        for run_str, lumis in payload.iteritems():
            run = int(run_str)
            for l in lumis:
                skimmer.AddGoodLumiRange(run, l[0], l[1])


# some common stuff that doesn't need to be configured
def run_Analyzer(skimmer, isData, output_name):
    # run and save output
    skimmer.Run()
    skimmer.Terminate()

    ret = path.isfile(output_name)
    if ret:
        logger.info(_sname+'.run_Analyzer', 'Successfully created %s'%(output_name))
        return output_name 
    else:
        logger.error(_sname+'.run_Analyzer', 'Failed in creating %s!'%(output_name))
        return False


def run_PandaAnalyzer(skimmer, isData, output_name):
    if isData:
        add_json(skimmer)

    return run_Analyzer(skimmer, isData, output_name)


def run_HRAnalyzer(*args, **kwargs):
    return run_Analyzer(*args, **kwargs) 

# main function to run a skimmer, customizable info 
# can be put in fn
def main(to_run, processed, fn):
    print_time('loading')
    for f in to_run.files:
        for i_attempt in xrange(maxcopy):
            input_name = request_data(f, i_attempt==0)
            if input_name is not None:
                break
            sleep(30)
        print_time('copy %s'%input_name)
        if input_name:
            logger.info(_sname+'.main',
                        'Starting to process '+input_name)
            success = fn(input_name, (to_run.dtype!='MC'), f)
            print_time('analyze %s'%input_name)
            if success:
                processed[input_name] = f
            if input_name[:5] == 'input': # if this is a local copy
                cleanup(input_name)
            print_time('remove %s'%input_name)
    
    if len(processed)==0:
        logger.warning(_sname+'.main', 'No successful outputs!')
        exit(1)



def wrapper(fn, pre_fn=None, post_fn=None):
    which = int(argv[1])
    submit_id = int(argv[2])

    sample_list = cb.read_sample_config('local.cfg',as_dict=False)
    to_run = None 
    for s in sample_list:
        if which==s.get_id():
            to_run = s
            break
    if not to_run:
        logger.error(_sname,'Could not find a job for PROCID=%i'%(which))
        exit(3)
    
    outdir = getenv('SUBMIT_OUTDIR')
    lockdir = getenv('SUBMIT_LOCKDIR')  
    outfilename = to_run.name+'_%i.root'%(submit_id)
    processed = {}
    
    report_start(outdir,outfilename,to_run.files)
    
    wd = isolate()
    if pre_fn is not None:
        pre_fn()
        print_time('pre_fn')

    main(to_run, processed, fn)
    
    hadd(processed.keys())
    print_time('hadd')

    if post_fn is not None:
        post_fn()
        print_time('post_fn')

    ret = stageout(outdir,outfilename)
    cleanup('*.root')
    un_isolate(wd)
    print_time('stageout and cleanup')
    if not ret:
        report_done(lockdir,outfilename,processed)
        cleanup('*.lock')
        print_time('create lock')
    else:
        exit(-1*ret)

    exit(0)
