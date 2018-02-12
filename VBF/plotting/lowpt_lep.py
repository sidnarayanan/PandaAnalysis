#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange
import numpy as np

basedir = getenv('PANDA_FLATDIR') + ''
infile = basedir+'/WJets.root'


parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cut',metavar='cut',type=str,default=None)
args = parser.parse_args()

figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Tools.Load import Load
from PandaCore.Tools.Misc import *
from PandaCore.Tools.root_interface import Selector
Load('HistogramDrawer')

lumi=35900
import PandaAnalysis.VBF.PandaSelection as sel 
base_cut = sel.cuts['signal']
label = 'inc'
if args.cut == 'cnc':
  base_cut = tAND(sel.cnc, sel.cuts['signal'])
  label = 'cnc'
elif args.cut == 'mjj':
  base_cut = tAND(sel.mjj, sel.cuts['signal'])
  label = 'mjj'
base_cut = 'pfmet>250  && dphipfmet>0.5 && nLooseMuon==0 && nLooseElectron==0 && nLoosePhoton==0 && fabs(calomet-pfmet)/pfmet<0.5 && nJet>0 && jet1Pt>100 && nTau==0 && jetNMBtags==0'
weight = sel.weights['signal']%lumi

plot = root.HistogramDrawer()
plot.SetRatioStyle()
plot.SetLumi(lumi/1000.)
plot.DrawEmpty(True)
plot.AddCMSLabel()
plot.AddLumiLabel()
#plot.SetAutoRange(False)
plot.InitLegend()
root.gStyle.SetOptStat(0)


s = Selector()
f = root.TFile.Open(infile); t = f.Get('events')
branches = ['pfmet', 'jot12Mass', '1',
            'fabs(genTauEta)', 'fabs(genMuonEta)', 'fabs(genElectronEta)', 
            'genTauPt', 'genElectronPt', 'genMuonPt',
            weight]
s.read_tree(t, branches = branches, cut = base_cut)

t_mask = ((s['genTauPt'] > 1) & (s['fabs(genTauEta)']<2.3) )
accp = ((s['genTauPt'] > 18))
t_out_mask = t_mask & (~accp) 


e_mask = ((s['genTauPt'] < 1) &
          (s['genElectronPt'] > 1) & (s['fabs(genElectronEta)']<2.4))
accp = ((s['genElectronPt'] > 10))
e_out_mask = e_mask & (~accp) 

m_mask = ((s['genTauPt'] < 1) &
          (s['genElectronPt'] < 1) &
          (s['genMuonPt'] > 1) & (s['fabs(genMuonEta)']<2.4))
accp = ((s['genMuonPt'] > 10))
m_out_mask = m_mask & (~accp) 


ratios_to_plot = {
  'pfmet' : ([250, 350, 500, 700, 1000],'E_{T}^{miss} [GeV]'),
  '1' : ([0, 2],'yield'),

  'jot12Mass' : ([200, 600, 1000, 1500, 2500, 3500, 5000],'m_{jj} [GeV]')
}
for k,v in ratios_to_plot.iteritems():
    h_inc = s.draw(k, weight, vbins=v[0])

    hs = {}

    hs['h_t_out'] = s.draw(k, weight, mask=t_out_mask, vbins=v[0])
    hs['h_m_out'] = s.draw(k, weight, mask=m_out_mask, vbins=v[0])
    hs['h_e_out'] = s.draw(k, weight, mask=e_out_mask, vbins=v[0])

    plot.Reset(False)
    for name,h in hs.iteritems(): 
        h.Divide(h_inc)
        h.GetXaxis().SetTitle(v[1])
        h.GetYaxis().SetTitle('Fraction of W events')
        h.SetMaximum(0.4)
        h.SetLineWidth(2)
        if k=='1':
          print name, h.Integral()

    plot.AddHistogram(hs['h_e_out'], 'W#rightarrowe p_{T} < 10 GeV', 7, -1, 'histe')
    plot.AddHistogram(hs['h_m_out'], 'W#rightarrow#mu p_{T} < 10 GeV', 6, -1, 'histe')
    plot.AddHistogram(hs['h_t_out'], 'W#rightarrow#tau p_{T} < 18 GeV', 8, -1, 'histe')


    plot.Draw(args.outdir,'ratio_'+label+'_'+k)
