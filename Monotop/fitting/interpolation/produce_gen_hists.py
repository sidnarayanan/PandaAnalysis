#!/usr/bin/env python

from glob import glob
from re import sub
from sys import argv,exit
from os import environ,system,path
from array import array

sname = argv[0].split('/')[-1]
m_V = int(argv[1])
m_DM = int(argv[2])
argv=[]

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import Load
from PandaCore.Tools.root_interface import read_tree, draw_hist, rename_dtypes

Load('Normalizer')

'''
    This script does the following things:
        - Read an input list of T2 files and xrdcp them locally
        - hadd these files
        - Remove the input files and only keep the merged one
        - Normalize the merged file
        - For each w in [1]+weights, draw the mediator pT distribution with weight (w*normalizedWeight)
        - Delete the merged file and stageout the file with histograms
'''

#CONFIG = 'couplings' # or finer or normal
CONFIG = 'normal'

list_dir = '/home/bmaier/cms/MonoTop/interpolation/'
xsec_path = 'non-resonant'
outdir = 'hists'
if CONFIG == 'couplings':
    list_dir += 'coupling_'
    xsec_path = 'non-resonant-couplings'
    outdir = 'hists_couplings'
elif CONFIG == 'fine':
    list_dir += 'fine_'
    xsec_path = 'non-resonant-fine'
    outdir = 'hists_fine'


## first copy files locally
def stage_in_file(source,target):
    source = 'root://xrootd.cmsaf.mit.edu/' + source
    source = source.replace('/mnt/hadoop/cms','')
    cmd = 'xrdcp %s %s'%(source,target)
    PInfo(sname+'.stage_in_file', cmd)
    system(cmd)

# copy slowly to keep Max happy
event_cutoff = 5e5
def stage_in_list():
    system('mkdir -p unmerged')
    flist = open(list_dir+'%i_%i.txt'%(m_V,m_DM))
    PInfo(sname+'.stage_in_list','Reading '+list_dir+'%i_%i.txt'%(m_V,m_DM))
    total_events = 0 
    for l in flist:
        in_name = l.strip()
        out_name = 'unmerged/'+in_name.split('/')[-1]
        stage_in_file(in_name,out_name)
        f_out = root.TFile(out_name)
        t_out = f_out.Get('events')
        total_events += t_out.GetEntries()
        if total_events >= event_cutoff:
            break

# do a recursive copy and make Max angry
def stage_in_files():
    system('mkdir -p unmerged')
    t2dir = open(list_dir+'%i_%i.txt'%(m_V,m_DM)).readlines()[0]
    t2dir = t2dir.split('/')[:-1]
    t2dir = '/'.join(t2dir)
    t2dir = t2dir.replace('/mnt/hadoop/cms','root://xrootd.cmsaf.mit.edu/')
    PInfo(sname+'.stage_in_files','Reading files from %s'%t2dir)
    cmd = 'xrdcp -r --parallel 4 %s unmerged'%(t2dir)
    system(cmd)

## now hadd and cleanup the files
def hadd():
    cmd = 'hadd -k -ff -n 100 -f merged.root unmerged/*root'
    system(cmd)

def remove(pattern):
    cmd = 'rm -rf %s'%pattern
    system(cmd)

## normalize the merged file
def get_xsec():
    params = read_nr_model(m_V,m_DM, path=xsec_path)
    if params:
        xsec = params.sigma
    else:
        exit(1)
    return xsec    

def normalize(xsec):
    f = root.TFile.Open('merged.root','UPDATE')
    t = f.Get('events')
    h = root.TH1D('h','',1,-1,2)
    t.Draw('1>>h','weight')
    total = h.Integral()

    n = root.Normalizer()
    n.isFloat = False
    n.inWeightName = 'weight'
    n.NormalizeTree(t,total,xsec)

    f.WriteTObject(t)
    f.Close()

## draw histograms
fweights = open(getenv('CMSSW_BASE')+'/src/PandaAnalysis/Monotop/fitting/signal_weights_all.dat')
if CONFIG == 'couplings':
    fweights = open(getenv('CMSSW_BASE')+'/src/PandaAnalysis/Monotop/fitting/signal_weights_couplings_all.dat')

weights = [x.strip() for x in fweights]
fweights.close()
bins = array('f', [175, 225, 275, 325, 375, 425, 475, 600, 800, 1200])
hbase = root.TH1D('hpt','',20,175,1200)
def draw_all():
    f_in = root.TFile.Open('merged.root')
    t_in = f_in.Get('events')
    f_out = root.TFile.Open('hists.root','RECREATE')
    weight_strs = []
    for idx in xrange(len(weights)):
        if idx==0:
            weight_str = 'normalizedWeight'
        else:
            weight_str = 'normalizedWeight*fabs(weights[%i])'%(idx-1)
        weight_strs.append(weight_str)

    n_weights = len(weight_strs)
    n_per = 50
    nominal_arr = None
    for iw in xrange(0, n_weights, n_per): 
        PInfo(sname,'Extracting %i -> %i'%(iw,iw+n_per))

        xarr = read_tree(t_in, ['genBosonPt'] + weight_strs[iw:iw + n_per])

        for idx in xrange(iw, iw + n_per):
            if idx == n_weights:
                break
            h = hbase.Clone()
            weight_str = weight_strs[idx]
            weight_name = weights[idx]
            draw_hist(h, xarr, ['genBosonPt'], weight_str)
            f_out.WriteTObject(h,'h_'+weight_name)

    f_out.Close()
    f_in.Close()

## stage out
def stageout():
    cmd = 'mkdir -p %s/interpolate/%s'%(getenv('PANDA_FITTING'),outdir)
    system(cmd)
    cmd = 'mv hists.root %s/interpolate/%s/%i_%i.root'%(getenv('PANDA_FITTING'),outdir,m_V,m_DM)
    PInfo(sname+'.stageout',cmd)
    system(cmd)

xsec = get_xsec() # do this first in case it's missing
#stage_in_files()
stage_in_list()
hadd()
remove('unmerged')
normalize(xsec)
draw_all()
remove('merged.root')
stageout()
remove('*root')
