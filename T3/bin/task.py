#!/usr/bin/env python

import json
import curses
import cPickle as pickle

import os 
from re import sub
from glob import glob
from query import query
from requests import post
from urllib2 import urlopen
from itertools import chain 
from time import time, sleep, strftime
from os import getenv,path,popen,makedirs
import subprocess as sp 
from shutil import rmtree
from itertools import chain

from PandaCore.Tools.script import * 
import PandaCore.Tools.job_management as jm

logger = Logger('task.py')

### Global definitions ###

validate_env()

# environment:
logdir        = getenv('SUBMIT_LOGDIR')
workdir       = getenv('SUBMIT_WORKDIR')
lockdir       = getenv('SUBMIT_LOCKDIR')
outdir        = getenv('SUBMIT_OUTDIR')
submit_name   = getenv('SUBMIT_NAME') + '_' + getenv('SUBMIT_USER')
cmssw_base    = getenv('CMSSW_BASE')
incfg         = workdir+'/local_all.cfg'
outcfg        = workdir+'/local.cfg'
panda_cfg     = getenv('PANDA_CFG')

args = parse(('--nfiles', {'type':int, 'default':-1}),
             ('--monitor', {'type':int, 'default':0}),
             *[(x, STORE_TRUE) for x in 
                ['--kill', '--kill_idle', '--check', '--submit', '--build_only', 
                 '--submit_only', '--clean_output', '--check_duplicates', '--rebuild', 
                 '--clean_duplicates', '--clean', '--force', '--silent', '--randomize']])

if args.clean:
    args.clean_output = True
if args.submit:
    args.submit_only = True
    args.build_only = True
if args.monitor:
    jm.SILENT = True
    args.check = True
if args.clean_duplicates:
    args.check_duplicates = True 
if args.check_duplicates and jm.textlock:
    logger.error('Duplicate checking is not yet implemented for text-locking, sorry!')
    sys.exit(1)

# for printing to screen:
columns = int(popen('stty size', 'r').read().split()[-1])
rows = int(popen('stty size', 'r').read().split()[0])
WIDTH = (columns-90)/2
header = ('%%-%is'%(WIDTH))%('Sample')
header += ('%%-%is'%(WIDTH+2))%('Progress')
header += ' %10s %10s %10s %10s %10s %10s'%('Total','T3', 'T2','Idle','Missing','Done')

colors = {
        'green' : 42,
        'blue' : 44,
        'cyan' : 46,
        'grey' : 47, 
        'red' : 41,
        }

last_lock = 1
last_check = 1

jm.setup_schedd(getenv('SUBMIT_CONFIG'))

def init_colors():
    curses.start_color()
    curses.init_pair(colors['green'], curses.COLOR_WHITE, curses.COLOR_GREEN)
    curses.init_pair(colors['blue'], curses.COLOR_WHITE, curses.COLOR_BLUE)
    curses.init_pair(colors['cyan'], curses.COLOR_WHITE, curses.COLOR_CYAN)
    curses.init_pair(colors['grey'], curses.COLOR_WHITE, curses.COLOR_WHITE)
    curses.init_pair(colors['red'], curses.COLOR_WHITE, curses.COLOR_RED)

### Task-specific functions ###

# for submission:
def build_snapshot(N, modify=False):
    if not modify:
        logger.info('Cleaning up staging areas...')
        rmtree(workdir); makedirs(workdir)
        rmtree(logdir); makedirs(logdir)
    else:
        logger.warning("Partially rebuilding work area on request...hope you know what you're doing")

    logger.info('Acquiring configuration...')
    do('wget -nv -O %s/list.cfg %s'%(workdir, panda_cfg))
    fin = open(outcfg.replace('local.cfg', 'list.cfg'))
    samples = jm.convert_catalog(list(fin), as_dict=True)
    keys = sorted(samples)
    if len(keys) == 0:
        logger.error('Configuration is empty - exiting!')
        sys.exit(1)

    to_write = [samples[k].get_config(N, suffix='_%i') for k in keys]
    with open(outcfg, 'w') as fout:
        for i,c in enumerate(chain.from_iterable(to_write)):
            fout.write(c%(i, i))
    logger.info('Submission will have %i jobs'%i)
    do('cp -v {0}/list.cfg {0}/list_all.cfg'.format(workdir))
    do('cp -v {0}/local.cfg {0}/local_all.cfg'.format(workdir))

    logger.info('Tarring up CMSSW...')
    cwd = os.getcwd()
    os.chdir(cmssw_base)
    do('tar --exclude-vcs -chzf cmssw.tgz src python biglib bin lib objs test external')
    do('mv -v cmssw.tgz %s'%workdir)

    logger.info('Creating executable...')
    os.chdir('%s/src/PandaAnalysis/T3/inputs/'%(cmssw_base))
    do('cp -v %s %s/skim.py'%(getenv('SUBMIT_TMPL'), workdir))
    do('chmod 775 %s/skim.py'%(workdir))

    logger.info('Finalizing work area...')
    jm.issue_proxy()
    do('cp -v /tmp/x509up_u%i %s/x509up'%(os.getuid(), workdir))
    do('cp -v %s/src/PandaAnalysis/T3/inputs/exec.sh %s'%(cmssw_base, workdir))

    logger.info('Taking a snapshot of work area...')
    do('cp -rvT %s %s/workdir'%(workdir, outdir))

    # input files for submission: cmssw.tgz, skim.py, x509up, local.cfg. exec.sh is the executable

def submit(silent=False):
    now = int(time())
    frozen_outcfg = outcfg.replace('local','local_%i'%now)
    do('cp %s %s'%(outcfg,frozen_outcfg)) 

    s = jm.Submission(frozen_outcfg,workdir+'/submission.pkl')
    s.execute(randomize=args.randomize)
    s.save()

    if not silent:
        statii = s.query_status()
        print 'Job summary:'
        for k,v in statii.iteritems():
            print '\t %10s : %5i'%(k,len(v))

def kill(idle=False):
    if path.isfile(workdir+'/submission.pkl'): 
        with open(workdir+'/submission.pkl','rb') as fpkl:
            submissions = pickle.load(fpkl)
            for s in submissions:
                s.kill(idle_only=idle)
            return
    else:
        logger.warning('Trying to kill a task with no submissions!')

def check_duplicates():
    url = jm.report_server + '/condor/query?task=%s'%(submit_name)
    completed = {}
    for r in json.load(urlopen(url)):
        completed.setdefault(r[0], []).append(r[1])
    to_clean = []
    for f,jids in completed.iteritems():
        if len(jids) > 1:
            logger.warning('%s found in %s'%(f, jids))
            to_clean += jids 
    if not to_clean:
        logger.info('No duplicates found')
    if args.clean_duplicates:
        for jid in to_clean:
            for f in glob('%s/*%s%'%(outdir, jid)):
                do('rm -f %s'%f)
            payload = {'task' : submit_name, 'job_id' : jid}
            post(jm.report_server+'/condor/clean', json=payload)


# for monitoring:
class Output:
    def __init__(self,name):
        self.name = name
        self.total = 0
        self.done = 0
        self.idle = 0
        self.t2 = 0
        self.t3 = 0
        self.missing = 0
    def add(self,state):
        self.total += 1
        setattr(self, state, getattr(self,state) + 1)
    def str(self):
        if self.total==0:
            return ''
        msgs = []
        s = ('%%-%is '%(WIDTH-1))%self.name[:(WIDTH-1)]
        if not args.monitor:
            msgs.append( s )
        else:
            msgs.append( (s,) )
        d_frac  = 1.*WIDTH*(self.done)/self.total
        t3_frac = 1.*WIDTH*(self.done+self.t3)/self.total
        t2_frac = 1.*WIDTH*(self.done+self.t3+self.t2)/self.total
        i_frac  = 1.*WIDTH*(self.done+self.t3+self.t2+self.idle)/self.total
        if not args.monitor:
            msgs.append( '[\033[0;%im'%colors['green'] )
        else:
            msgs.append( '[' )
        state = 0
        for i in xrange(WIDTH):
            if not args.monitor:
                s = ''
                if i>=i_frac:
                    if state != 4:
                        s += '\033[0;%im'%colors['red']
                        state = 4
                elif i>=t2_frac:
                    if state!=3:
                        s += '\033[0;%im'%colors['grey']
                        state = 3
                elif i>=t3_frac:
                    if state!=2:
                        s += '\033[0;%im'%colors['cyan']
                        state = 2
                elif i>=d_frac:
                    if state!=1:
                        s += '\033[0;%im'%colors['blue']
                        state = 1
                s += ' '
                msgs.append( s )
            else:
                if i>=i_frac:
                    msgs.append( (' ', curses.color_pair(colors['red'])) )
                elif i>=t2_frac:
                    msgs.append( (' ', curses.color_pair(colors['grey'])) )
                elif i>=t3_frac:
                    msgs.append( (' ', curses.color_pair(colors['cyan'])) )
                elif i>=d_frac:
                    msgs.append( (' ', curses.color_pair(colors['blue'])) )
                else:
                    msgs.append( (' ', curses.color_pair(colors['green'])) )
        s = ''
        if not args.monitor:
            s += '\033[0m] '
        else:
            s += ']'
        s += '%10i '%self.total
        s += '%10i '%self.t3
        s += '%10i '%self.t2
        s += '%10i '%self.idle
        s += '%10i '%self.missing
        s += '%10i '%self.done
        s += '(done=%.2f%%)\n'%(d_frac*100./WIDTH)
        if not args.monitor:
            msgs.append( s )
            return ''.join(msgs)
        else:
            msgs.append( (s,) )
            return msgs

def check(stdscr=None):
    if args.monitor:
        init_colors()
        stdscr.nodelay(True)
    while True:
        global last_lock
        refresh = False
        if args.monitor:
            c = stdscr.getch()
            curses.flushinp()
            if c == ord('q'):
                return
            refresh = (c == ord('r'))

        if time() - args.monitor > last_lock:
            refresh = True 

        if not refresh:
            sleep(1)
            continue 

        processedfiles = set([])
        if jm.textlock:
            # determine what files have been processed and logged as such
            locks = glob(lockdir+'/*lock')
            nl = len(locks)
            il = 1
            for lock in locks:
                il += 1
                try:
                    flock = open(lock)
                    for l in flock:
                        processedfiles.add(l.strip())
                except IOError:
                    pass
        else:
            url = jm.report_server + '/condor/query?task=%s'%(submit_name)
            for r in json.load(urlopen(url)):
                if r[2] is not None and r[-1] == 0:
                    processedfiles.add(r[0])

        # determine what samples from previous resubmissions are still running
        t2_samples = []
        t3_samples = []
        idle_samples = []
        if path.isfile(workdir+'/submission.pkl'): 
            with open(workdir+'/submission.pkl','rb') as fpkl:
                submissions = pickle.load(fpkl)
        else:
            submissions = []
        for s in submissions:
            results = s.query_status()
            t3_samples += results['T3']
            t2_samples += results['T2']
            idle_samples += results['idle']

        t2_files = list(chain.from_iterable([x.files for x in t2_samples]))
        t3_files = list(chain.from_iterable([x.files for x in t3_samples]))
        idle_files = list(chain.from_iterable([x.files for x in idle_samples]))

        # for fancy display
        outputs = {}
        data = Output('Data')
        mc = Output('MC')

        all_samples = jm.read_sample_config(incfg)
        filtered_samples = {}
        merged_samples = {}
        outfile = open(outcfg,'w')
        for name in sorted(all_samples):
            sample = all_samples[name]
            out_sample = jm.DataSample(name,sample.dtype,sample.xsec)

            base_name = sub('_[0-9]+$','',name)
            if base_name not in outputs:
                outputs[base_name] = Output(base_name)
            output = outputs[base_name]
            if base_name not in merged_samples:
                merged_samples[base_name] = jm.DataSample(base_name,sample.dtype,sample.xsec)
            merged_sample = merged_samples[base_name]

            to_resubmit = []

            for f in sample.files:
                state = 'missing'
                hashed = f if jm.textlock else jm.md5hash(f) 
                if hashed in processedfiles:
                    state = 'done'
                elif f in t3_files:
                    state = 't3'
                elif f in t2_files:
                    state = 't2'
                elif f in idle_files:
                    state = 'idle'

                if state=='missing' or (args.force and state!='done'):
                    out_sample.add_file(f)
                    merged_sample.add_file(f)

                output.add(state)
                if sample.dtype=='MC':
                    mc.add(state)
                else:
                    data.add(state)

            if len(out_sample.files)>0:
                filtered_samples[name] = out_sample

        if args.nfiles<0:
            keys = sorted(filtered_samples)
            for k in keys:
                sample = filtered_samples[k]
                if len(sample.files)==0:
                    continue
                configs = sample.get_config(-1)
                for c in configs:
                    outfile.write(c)
        else:
            keys = sorted(merged_samples)
            counter=0
            for k in keys:
                sample = merged_samples[k]
                if len(sample.files)==0:
                    continue
                configs = sample.get_config(args.nfiles,suffix='_%i')
                for c in configs:
                    outfile.write(c%(counter,counter))
                    counter += 1

        msg = ['TASK = '+submit_name]
        if args.monitor:
            msg.append('\n')

        msg.append(header)
        if args.monitor:
            msg.append('\n')

        if args.monitor and len(outputs)+10>rows:
            args.silent = True
            msg.append( ('Too many samples to show in monitoring mode!\n', curses.color_pair(colors['red'])) )
        if not args.silent:
            for n in sorted(outputs):
                if args.monitor:
                    msg.extend(outputs[n].str())
                else:
                    msg.append(outputs[n].str().strip())
            msg.append('')
        if args.monitor:
            msg.extend(data.str())
            msg.extend(mc.str())
        else:
            msg.append(data.str().strip())
            msg.append(mc.str().strip())
        msg.append('')
        if args.monitor:
            msg.append('Legend: Done=[')
            msg.append( ('    ',curses.color_pair(colors['green'])) )
            msg.append('], T3=[')
            msg.append( ('    ',curses.color_pair(colors['blue'])) )
            msg.append('], T2=[')
            msg.append( ('    ',curses.color_pair(colors['cyan'])) )
            msg.append('], Idle[')
            msg.append( ('    ',curses.color_pair(colors['grey'])) )
            msg.append('], Missing=[')
            msg.append( ('    ',curses.color_pair(colors['red'])) )
            msg.append(']\n')
        else:
            msg.append( 'Legend: Done=\033[0;%im    \033[0m, T3=\033[0;%im    \033[0m, T2=\033[0;%im    \033[0m, Idle=\033[0;%im    \033[0m, Missing=\033[0;%im    \033[0m, '%(colors['green'],colors['blue'],colors['cyan'],colors['grey'],colors['red']))

        outfile.close()

        msg.append(strftime('%Y-%m-%d %H:%M:%S'))
        if args.monitor:
            msg.append('\n')
        msg.append( '\nMost recent submission:')
        if args.monitor:
            msg.extend([x+'\n' for x in query()])
            msg.append('\nPress "r" to refresh or "q" to close')
        else:
            msg.extend(query())
        msg.append('')

        if args.monitor:
            stdscr.clear()
            for m in msg:
                if type(m) == str:
                    stdscr.addstr(m)
                else:
                    stdscr.addstr(*m)
            stdscr.refresh()
        else:
            sys.stdout.write('\n'.join(msg))

        last_lock = int(time())

        if args.submit_only and (mc.missing + data.missing > 0):
            submit(silent=(args.monitor is not None))


        if args.monitor:
            sleep(1)
        else:
            return

### MAIN ###
if not args.monitor:
    logger.info('TASK = '+submit_name)

if args.check_duplicates:
    check_duplicates()

if args.kill or args.kill_idle:
    kill(args.kill_idle)

if args.clean_output:
    if jm.textlock:
        logger.info('Cleaning up %s and %s'%(lockdir, outdir))
        sleep(2)
        do('rm -rf %s %s'%(lockdir, outdir))
        do('mkdir -p %s %s'%(lockdir, outdir))
    else:
        logger.info('Cleaning up %s and deleting entries'%(outdir))
        sleep(2)
        do('rm -rf %s'%(outdir))
        do('mkdir -p %s'%(outdir))
        payload = {'task' : submit_name}
        post(jm.report_server+'/condor/clean', json=payload)
    if args.clean:
        logger.info('Cleaning up %s and %s'%(logdir, workdir))
        sleep(2)
        do('rm -rf %s %s'%(logdir, workdir))
        do('mkdir -p %s %s'%(logdir, workdir))

if args.check:
    if args.monitor:
        curses.wrapper(check)
    else:
        check()
else:
    if (args.build_only or args.rebuild) and (not path.isfile(workdir+'/submission.pkl') or not args.submit): 
        if args.nfiles < 0:
            logger.info('Number of files not provided for new task => setting nfiles=25')
            args.nfiles = 25
        build_snapshot(args.nfiles, modify=args.rebuild)
    if args.submit_only:
        submit(silent=False)
