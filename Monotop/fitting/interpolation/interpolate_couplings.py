#!/usr/bin/env python
from re import sub
from sys import argv,exit
from os import path,getenv,system
from array import array
from glob import glob
import argparse

parser = argparse.ArgumentParser(description='interpolate')
parser.add_argument('--mass',type=str)
parser.add_argument('--mass_reco',type=str,default=None)
args = parser.parse_args()
sname = argv[0].split('/')[-1]
argv=[]

from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import Load
from PandaCore.Tools.root_interface import read_tree, draw_hist
import ROOT as root
from math import sqrt

Load('BranchAdder')
ba = root.BranchAdder()
ba.verbose = False
ba.formula = 'genBosonPt'

bins = array('f',[250,280,310,350,400,450,600,1000])
N = len(bins)-1
h_base = root.TH1F('hbase','',N,bins)

systs = ['btagUp','btagDown','mistagUp','mistagDown',
         'sjbtagUp','sjbtagDown','sjmistagUp','sjmistagDown']


basedir = getenv('PANDA_FITTING')

couplings_interp = map(lambda x : x.strip(), 
                       open(getenv('CMSSW_BASE')+'/src/PandaAnalysis/Monotop/fitting/signal_weights_couplings_all.dat').readlines())
# couplings_reco = map(lambda x : x.strip(), 
#                      open(getenv('CMSSW_BASE')+'/src/PandaAnalysis/Monotop/fitting/signal_weights_all.dat').readlines())

def parse_coupling_name(x):
    if x=='nominal':
        return (1, 0, 0.25, 0)
    x = x.replace('rw_','').replace('_nlo','')
    xx = x.split('_')
    gdmv = float(xx[1].replace('p','.'))
    gdma = float(xx[3].replace('p','.'))
    gv = float(xx[5].replace('p','.'))
    ga = float(xx[7].replace('p','.'))
    return (gdmv, gdma, gv, ga)

def euclid(x,y):
    dim = len(x)
    r = 0
    for i in xrange(dim):
        r += pow(x[i] - y[i], 2)
    r = sqrt(r)
    return r

# todo: pass more complicated metrics to interpolate, including things
# that can take N nearby points and weight them
def interpolate():
    f_interp_path = basedir+'/interpolate/hists_couplings/%s.root'%(args.mass)
    f_interp = root.TFile.Open(f_interp_path)

    PInfo(sname, 'We want to interpolate point %s'%args.mass)

    if not args.mass_reco:
        m_V,m_DM = [int(x) for x in args.mass.split('_')]
        offshell = (2*m_DM >= m_V)
        frescos = map(lambda x : x.split('/')[-1].replace('.root','').replace('fittingForest_signal_vector_',''),
                      glob(basedir+'/signals_3/fittingForest_signal_vector*root'))
        frescos = map(lambda x : tuple(map(lambda y : int(y), x.split('_'))),
                      frescos)
        m_V_recos = [x[0] for x in frescos]
        m_V_reco = min(m_V_recos, key = lambda x : abs(m_V-x))

        m_DM_recos = [x[1] for x in frescos if x[0]==m_V_reco]
        m_DM_recos = [x for x in m_DM_recos if ((2*x>=m_V)==offshell)]
        if len(m_DM_recos)==0:
            m_DM_recos = [x[1] for x in frescos if x[0]==m_V_reco]
        m_DM_reco = min(m_DM_recos, key = lambda x : abs(m_DM-x))

        args.mass_reco = '%i_%i'%(m_V_reco,m_DM_reco)

    PInfo(sname, 'We have chosen as the reconstructed point %s'%args.mass_reco)

    f_reco_path = basedir+'/signals_3/fittingForest_signal_vector_%s.root'%(args.mass_reco)
    f_reco = root.TFile.Open(f_reco_path)
    t_reco = f_reco.Get('%s_signal'%(args.mass_reco))
    f_reco_gen_path = basedir+'/interpolate/hists/%s.root'%(args.mass_reco)
    f_reco_gen = root.TFile.Open(f_reco_gen_path)

    system('mkdir -p %s/interpolate/interpolated_couplings2test/'%basedir)
    t_out = t_reco.Clone()
    metname = 'min(met,999.9999)'
    branches = [metname]
    interp2reco = {}
    for g in couplings_interp:
        # find the reco g
        # g_reco = min(couplings_reco, key = lambda x : euclid(parse_coupling_name(g), parse_coupling_name(x)))
        g_reco = 'nominal' # override
        PInfo(sname,'Interpolating %s with %s'%(g,g_reco))
        interp2reco[g] = g_reco
        h_interp = f_interp.Get('h_%s'%g)
        h_reco_gen = f_reco_gen.Get('h_%s'%g_reco)
        h_ratio = h_interp.Clone(); h_ratio.Divide(h_reco_gen)
        ba.newBranchName = 'interp_%s'%g
        ba.AddBranchFromHistogram(t_reco,h_ratio)
        if g_reco=='nominal':
            weightname = 'weight*interp_%s'%g
        else:
            weightname = 'weight*%s*interp_%s'%(g_reco.replace('_nlo','').replace('rw_',''),g)
        branches.append(weightname)
        for s in systs:
            branches.append(weightname.replace('weight',s))
    xarr = {}
    xarr['tight'] = read_tree(t_reco,branches,'top_ecf_bdt>0.45 && met>250')
    xarr['loose'] = read_tree(t_reco,branches,'top_ecf_bdt>0.1 && top_ecf_bdt<0.45 && met>250')

    f_out = root.TFile('%s/interpolate/interpolated_couplings2test/interp_%s.root'%(basedir,args.mass),'RECREATE')
    PInfo(sname,'Creating %s/interpolate/interpolated_couplings2test/interp_%s.root'%(basedir,args.mass))
    for g in couplings_interp:
        f_out.cd()
        folder = f_out.mkdir(g)
        h_templates = {}
        g_reco = interp2reco[g]
        if g_reco=='nominal':
            weightname = 'weight*interp_%s'%g
        else:
            weightname = 'weight*%s*interp_%s'%(g_reco.replace('_nlo','').replace('rw_',''),g)
        for t in ['tight','loose']:
            h_templates[t] = h_base.Clone()
            draw_hist(h_templates[t], xarr[t], [metname], weightname)
            for s in systs:
                hname = '%s_%s'%(t,s)
                h_templates[hname] = h_base.Clone()
                draw_hist(h_templates[hname], xarr[t], [metname], 
                          weightname.replace('weight',s))
        for k,h in h_templates.iteritems():
            folder.WriteTObject(h,'h_%s'%(k))
    f_out.Close()


interpolate()
