#!/usr/bin/env python

from sys import argv, exit
import numpy as np 
from os import getenv, system
from PandaCore.Tools.Misc import PInfo, PDebug 
import PandaAnalysis.Deep.job_deepgen_utilities as deep_utils
from glob import glob 

deep_utils.NORM = True
truth = ['nprongs']
fractions = {'train':0.6, 'test':0.15}
fcfg = open(argv[1])
name = argv[2]
singletons = None
outdir = getenv('SUBMIT_NPY')
datadir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/deep/'
submit_name = getenv('SUBMIT_NAME')
me = argv[0].split('/')[-1]
argv = []

n_partons_proc = {
            'QCD'   : 1,
            'Top'   : 3,
            'Top_lo': 3,
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
f_pt = root.TFile.Open(datadir + 'flatten_gen.root')
h_pt = f_pt.Get('h_%s'%('_'.join(name.split('_')[:-1])))
f_pt_scaled = root.TFile.Open(datadir + 'flatten_gen_scaled.root')
h_pt_scaled = f_pt_scaled.Get('h_%s'%('_'.join(name.split('_')[:-1])))

print 'h_%s_%s'%('_'.join(name.split('_')[:-1]), submit_name)

data = {}
for fpath in fcfg.readlines():
    print fpath.strip()
    try:
        d = np.load(fpath.strip())
    except [IOError, AttributeError] as e:
        PError(me, str(e))
        continue
    mask = (d['nprongs'] == n_partons)
    for k,v in d.iteritems():
        if k == 'singleton_branches':
            data[k] = v 
            continue
        if v.shape[0]:
            if k not in data:
                data[k] = []
            #data[k].append(v)
            data[k].append(v[mask])

if not len(data):
    PInfo(me, 'This was an empty config!')
    exit(0)

for k,v in data.iteritems():
    if k == 'singleton_branches':
        continue
    data[k] = np.concatenate(v)

if not data['pt'].shape[0]:
    PInfo(me, 'Nothing passed the mask')
    exit(0)

if deep_utils.NORM:
    deep_utils.normalize_arrays(data, 'particles')


def reweight(x_pt):
    return h_pt.GetBinContent(h_pt.FindBin(x_pt))
reweight = np.vectorize(reweight)

def reweight_s(x_pt):
    return h_pt_scaled.GetBinContent(h_pt_scaled.FindBin(x_pt))
reweight_s = np.vectorize(reweight_s)


data['ptweight'] = reweight(data['pt'])
data['ptweight_scaled'] = reweight_s(data['pt'])

def xform(x, pt, name):
    if re.match(r'^[0-9]_[0-9]_[0-9]$', name):
        exponent = int(name.split('_')[0])
        return x / np.power(pt, exponent)
    else:
        return x

def dump(idx, partition):
    global singletons 

    outpath = 'tmp/' + partition + '/' + name + '_%s.npy'

    # singletons
    singletons = data['singleton_branches']
    d = np.vstack([xform(data[x][idx], data['pt'], x) for x in singletons]).T 
    np.save(outpath%'singletons', d)

    # particles
    d = data['particles'][idx, :, :]
    np.save(outpath%'particles', d)

    # truth
    d = np.vstack([data[x][idx] for x in truth]).T 
    np.save(outpath%'truth', d)

    # pt weights
    d = data['ptweight'][idx]
    np.save(outpath%'ptweight', d)

    d = data['ptweight_scaled'][idx]
    np.save(outpath%'ptweight_scaled', d)
    


pt = data['pt']
mask = np.logical_and(pt > 450, pt < 1200)

indices = np.array(range(data['pt'].shape[0]))
indices = indices[mask] # only within pT window
np.random.shuffle(indices)

N = {k:int(len(indices) * v) for k,v in fractions.iteritems()}

for d in ['train', 'test', 'validate']:
    system('mkdir -p tmp/%s'%d)

dump(indices[:N['train']], 'train')
dump(indices[N['train']:N['train']+N['test']], 'test')
dump(indices[N['train']+N['test']:], 'validate')

print "singletons = ["
for s in singletons:
    print "    '%s',"%s
print "]"

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
