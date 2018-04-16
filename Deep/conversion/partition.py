#!/usr/bin/env python

from glob import glob
from os import getenv, system, path
from argparse import ArgumentParser
from re import sub
from random import shuffle

parser = ArgumentParser()
parser.add_argument('--nmax',type=int,default=25)
parser.add_argument('--nmin',type=int,default=10)
parser.add_argument('--proc',type=str)
args = parser.parse_args()

outdir = getenv('SUBMIT_OUTDIR')
npydir = getenv('SUBMIT_NPY')

if args.proc == 'Top':
    fs = []
    for p in ['ZpTT', 'Scalar_MonoTop', 'Vector_MonoTop']:
        fs +=  glob(outdir + '/' + p + '*.npz')
elif args.proc == 'Top_lo':
    fs = []
    for p in ['ZpTT_lo']:
        fs +=  glob(outdir + '/' + p + '*.npz')
elif args.proc == 'Higgs':
    fs = []
    for p in ['ZpA0h']:
        fs +=  glob(outdir + '/' + p + '*.npz')
elif args.proc == 'W':
    fs = []
    for p in ['ZpWW']:
        fs +=  glob(outdir + '/' + p + '*.npz')
elif 'QCD_' in args.proc:
    fs = []
    fs = glob(outdir + '/' + args.proc.replace('QCD_','QCD_*_') + '*.npz')
else:
    fs = glob(outdir + '/' + args.proc + '*.npz')

pd_map = {}
for f in fs:
    f_ = f.split('/')[-1]
    pd = sub('_[0-9]+_[0-9]+.npz','',f_)
    if pd not in pd_map:
        pd_map[pd] = []
    pd_map[pd].append(f)


npartition = max(args.nmin, min(map(len, pd_map.values())))

pdir = npydir + '/partitions/'
system('mkdir -p ' + pdir + args.proc)
system('rm -rf ' + pdir + args.proc + '/*')

arglist = []

for k in xrange(npartition):
    to_run = []
    for pd,fs in pd_map.iteritems():
        n_this = max(1, len(fs) / npartition)
        lo = n_this * k
        hi = -1 if (k == npartition - 1) else n_this * (k + 1)
        to_run += fs[lo:hi]

    shuffle(to_run)

    if len(to_run) <= args.nmax:
        klabel = str(k)
        with open(pdir + args.proc + '/' + klabel + '.txt', 'w') as fout:
            for r in to_run:
                fout.write(r + '\n')
        arglist.append( pdir + args.proc + '/' + klabel + '.txt' )
    else:
        nsub = len(to_run) / args.nmax + 1
        for kk in xrange(nsub):
            lo = args.nmax * kk 
            hi = -1 if (kk == nsub - 1) else args.nmax * (kk + 1)
            klabel = str(k) + '_' + str(kk)
            with open(pdir + args.proc + '/' + klabel + '.txt', 'w') as fout:
                for r in to_run[lo:hi]:
                    fout.write(r + '\n')
            arglist.append( pdir + args.proc + '/' + klabel + '.txt' )

shuffle(arglist)
with open(pdir + args.proc + '.txt', 'w') as fout:
    for i,a in enumerate(arglist):
        fout.write(path.realpath(a) + ' %s_%i \n'%(args.proc,i))
print len(arglist)
