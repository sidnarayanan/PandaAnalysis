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
from PandaCore.Tools.script import * 

## TODO: matrix of errors correlated with
## hosts where job was running

logdir = getenv('SUBMIT_LOGDIR')
workdir = getenv('SUBMIT_WORKDIR')
args = parse(('--submission', {'type':int, 'nargs':'+', 'default':None}),
             ('--verbose', STORE_TRUE),
             ('--dump', STORE_TRUE),
             ('--nodone', STORE_TRUE))
logdirpath = '$SUBMIT_LOGDIR/' if args.submission is None else '$SUBMIT_LOGDIR/[%s]_'%(''.join(map(str,args.submission)))
dumpdirpath = getenv('SUBMIT_LOGDIR')+'/logs_dump/'

def sub(x, y, z):
    return rsub(y, z, x)

def clean(x):
    x = x.replace('+', '\+')
    x = sub(x, r'\n', '')
    x = sub(x, r'[\.0-9A-Za-z\-_]*\.root', 'X.root')
    x = sub(x, r'[0-9][0-9][0-9][0-9]+', 'X')
    x = sub(x, r'/mnt/hadoop.*root', 'X.root')
    x = sub(x, r'/store/user.*root', 'X.root')
    x = sub(x, r'/data/t3.*lock', 'X.lock')
    x = sub(x, r'branch:.*', '')
    x = sub(x, r'/mnt/hadoop.*npz', 'X.npz')
    x = sub(x, r'object at [a-zA-Z0-9]*', 'object at ADDRESS') 
    x = sub(x, r'\033\[91m', '')
    x = sub(x, r't3btch[0-9][0-9][0-9]', 't3btchX')
    return x 

cmd = 'grep -i "error\|fatal" %s*err'%logdirpath
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
    do('rm -rf %s'%dumpdirpath)
    do('mkdir -p %s'%dumpdirpath)

done = set([])
cmd = r'grep -i "report_done.*Response \[200\]"  %s*err'%logdirpath
for l in subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.readlines():
    f = sub(l.split(':')[0].split('/')[-1], '.err', '')
    done.add(f)

for i,c in enumerate(correlations):
    if args.verbose:
        print 'Failed with error class %i:'%i
    files = {}
    hosts = {}
    for f in sorted(aggregates[list(c)[0]]):
        if args.nodone and f in done:
            continue
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
                    hosts[ll] = set([])
                hosts[ll].add( '%s/%s.err'%(logdir, f) )
    if args.verbose:
        for f,n in files.iteritems():
            print '   ',f,n
    else:
        try:
            print 'Failure class %2i failed on %3i files, an average of %.1f times'%(i, len(files), float(sum(files.values())) / len(files) / 2) # each file appears in error logs twice
        except ZeroDivisionError:
            print 'Failure class %2i failed on %3i jobs, but number of files is unknown'%(i, len(aggregates[list(c)[0]])) 
    if args.dump:
        with open('%s/%i.log'%(dumpdirpath,i), 'w') as fdump:
            fdump.write('Error summary:\n')
            for msg in c:
                fdump.write('\t' + msg + '\n')
            for h,v in hosts.iteritems():
                fdump.write(h + ':\n')
                for vv in v:
                    fdump.write('\t' + vv + '\n')

