#!/usr/bin/env python

from os import getenv,path,popen,system
from PandaCore.Tools.job_management import *
import subprocess
import sys
import argparse
from glob import glob
from re import sub as rsub
import cPickle as pickle
from itertools import chain 

## TODO: matrix of errors correlated with
## hosts where job was running

logdir = getenv('SUBMIT_LOGDIR')
workdir = getenv('SUBMIT_WORKDIR')
parser = argparse.ArgumentParser(description='analyze log files')
parser.add_argument('--verbose',action='store_true')
parser.add_argument('--dump',action='store_true')
parser.add_argument('--submission',type=int,nargs='+',default=None)
args = parser.parse_args()
logdirpath = '$SUBMIT_LOGDIR/' if args.submission is None else '$SUBMIT_LOGDIR/[%s]_'%(''.join(map(str,args.submission)))

def sub(x, y, z):
    return rsub(y, z, x)

def clean(x):
    x = x.replace('+', '\+')
    x = sub(x, '\n', '')
    x = sub(x, 'input_[.0-9A-Z\-]*\.root', 'X.root')
    x = sub(x, '[0-9][0-9][0-9][0-9]+', 'X')
    x = sub(x, '/mnt/hadoop.*root', 'X.root')
    x = sub(x, '/store/user.*root', 'X.root')
    x = sub(x, '/data/t3.*lock', 'X.lock')
    x = sub(x, 'branch:.*', '')
    x = sub(x, '/mnt/hadoop.*npz', 'X.npz')
    return x 

cmd = 'grep -i error %s*err'%logdirpath
errors_by_file = {}
for l in subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.readlines():
    f = sub(l.split(':')[0].split('/')[-1], '.err', '')
    if f not in errors_by_file:
        errors_by_file[f] = []
    msg = clean(l.replace(l.split(':')[0], '')[1:])
    errors_by_file[f].append(msg)


aggregates = {}
for f,errors in errors_by_file.iteritems():
    for e in errors:
        if e not  in aggregates:
            aggregates[e] = set([f])
        else:
            aggregates[e].add(f)


correlations = []
sorted_aggregates = sorted(aggregates)
for i in xrange(len(sorted_aggregates)):
    a = sorted_aggregates[i]
    added = False
    for j in xrange(i):
        b = sorted_aggregates[j]
        if aggregates[a] == aggregates[b]:
            for c in correlations:
                if b in c:
                    c.add(a)
                    added = True
                    break
    if not added:
      correlations.append(set([a]))

cache = pickle.load(open(getenv('SUBMIT_WORKDIR') + '/submission.pkl', 'rb'))


for i,c in enumerate(correlations):
    print 'Failure class %i:'%i
    for msg in c:
        print '   ',msg

if args.dump:
    system('mkdir -p logs_dump')
    system('rm -f logs_dump/*')

for i,c in enumerate(correlations):
    if args.verbose:
        print 'Failed with error class %i:'%i
    files = {}
    hosts = {}
    for f in sorted(aggregates[list(c)[0]]):
        cmd = 'grep -o "pandaf[^ ]*\.root" $SUBMIT_LOGDIR/%s.err'%(f)
        for l in subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.readlines():
            ll = l.strip()
            if ll not in files:
                files[ll] = 0
            files[ll] += 1
        if args.dump:
            cmd = 'grep -o "hostname = .*" $SUBMIT_LOGDIR/%s.err'%(f)
            for l in subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.readlines():
                ll = l.strip()
                if ll not in hosts:
                    hosts[ll] = []
                hosts[ll].append( '%s/%s.err'%(logdir, f))
    if args.verbose:
        for f,n in files.iteritems():
            print '   ',f,n
    else:
        try:
            print 'Failure class %2i failed on %3i files, an average of %.1f times'%(i, len(files), float(sum(files.values())) / len(files) / 2) # each file appears in error logs twice
        except ZeroDivisionError:
            print 'Failure class %2i failed on %3i jobs, but number of files is unknown'%(i, len(aggregates[list(c)[0]])) 
    if args.dump:
        with open('logs_dump/class_%i.log'%i, 'w') as fdump:
            fdump.write('Error summary:\n')
            for msg in c:
                fdump.write('\t' + msg + '\n')
            for h,v in hosts.iteritems():
                fdump.write(h + ':\n')
                for vv in v:
                    fdump.write('\t' + vv + '\n')

