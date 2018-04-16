#!/usr/bin/env python

from sys import argv, exit
import numpy as np 
from os import getenv, system
from PandaCore.Tools.Misc import PInfo, PDebug 
import PandaAnalysis.Deep.job_deep_utilities as deep_utils
from glob import glob 

NORM = True
singletons = deep_utils.singleton_branches
truth = ['nPartons', 'nBPartons', 'nCPartons']
events = ['eventNumber']
fractions = {'train':0.7, 'test':0.15}
fcfg = open(argv[1])
name = argv[2]
outdir = getenv('SUBMIT_NPY')
datadir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/deep/'
me = argv[0].split('/')[-1]
argv = []

n_partons_proc = {
            'QCD'   : 1,
            'Top'   : 3,
            'ZpTT'  : 3,
            'ZpWW'  : 2,
            'ZpA0h' : 2,
            'Higgs' : 2,
            'W'     : 2,
        }
n_partons = 1
for k,v in n_partons_proc.iteritems():
    if k in name:
        n_partons = v
        break


import ROOT as root 
f_pt = root.TFile.Open(datadir + 'flatten.root')
h_pt = f_pt.Get('h_%s'%(name.split('_')[0]))
f_pt_scaled = root.TFile.Open(datadir + 'flatten_scaled.root')
h_pt_scaled = f_pt_scaled.Get('h_%s'%(name.split('_')[0]))

data = {}
for fpath in fcfg.readlines():
    d = np.load(fpath.strip())
    mask = (d['nPartons'] == n_partons)
    for k,v in d.iteritems():
        if v.shape[0]:
            if k not in data:
                data[k] = []
            data[k].append(v[mask])

if not len(data):
    PInfo(me, 'This was an empty config!')
    exit(0)

for k,v in data.iteritems():
    data[k] = np.concatenate(v)

if not data['pt'].shape[0]:
    PInfo(me, 'Nothing passed the mask')
    exit(0)

if NORM:
    deep_utils.normalize_arrays(data, 'pf')
    deep_utils.normalize_arrays(data, 'sv')


def reweight(x_pt):
#    x_pt = 400 + (600 * x_pt)
    return h_pt.GetBinContent(h_pt.FindBin(x_pt))
reweight = np.vectorize(reweight)

def reweight_s(x_pt):
#    x_pt = 400 + (600 * x_pt)
    return h_pt_scaled.GetBinContent(h_pt_scaled.FindBin(x_pt))
reweight_s = np.vectorize(reweight_s)


data['ptweight'] = reweight(data['rawpt'])
data['ptweight_scaled'] = reweight_s(data['rawpt'])


def dump(idx, partition):
    outpath = 'tmp/' + partition + '/' + name + '_%s.npy'

    # singletons
    d = np.vstack([data[x][idx] for x in singletons]).T 
    np.save(outpath%'singletons', d)

    # events
    d = np.vstack([data[x][idx] for x in events]).T 
    np.save(outpath%'events', d)

    # pf
    d = data['pf'][idx, :, :]
    np.save(outpath%'pf', d)

    # sv
    d = data['sv'][idx, :, :]
    np.save(outpath%'sv', d)

    # truth
    d = np.vstack([data[x][idx] for x in truth]).T 
    np.save(outpath%'truth', d)

    # pt weights
    d = data['ptweight'][idx]
    np.save(outpath%'ptweight', d)

    d = data['ptweight_scaled'][idx]
    np.save(outpath%'ptweight_scaled', d)
    


pt = data['rawpt']
mask = np.logical_and(pt > 450, pt < 1200)

indices = np.array(range(data['eventNumber'].shape[0]))
indices = indices[mask] # only within pT window
np.random.shuffle(indices)

N = {k:int(len(indices) * v) for k,v in fractions.iteritems()}

for d in ['train', 'test', 'validate']:
    system('mkdir -p tmp/%s'%d)

dump(indices[:N['train']], 'train')
dump(indices[N['train']:N['train']+N['test']], 'test')
dump(indices[N['train']+N['test']:], 'validate')

ret = None
for d in ['train', 'test', 'validate']:
    for ftmp in glob('tmp/'+d+'/*npy'):
        cmd = 'cp -v %s %s/%s'%(ftmp,outdir,ftmp.replace('tmp/',''))
        # cmd = 'gfal-copy -f file://$PWD/%s srm://t3serv006.mit.edu:8443/srm/v2/server?SFN=%s/%s'%(ftmp,outdir,ftmp.replace('tmp/',''))
        PInfo(me, cmd)
        ret = max(ret, system(cmd))

system('rm -rf tmp')
PDebug(me, 'exit code %i'%ret)
exit(ret)
