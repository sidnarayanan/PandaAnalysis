#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange
import numpy as np

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--pd',metavar='pd',type=str,default='MET')
args = parser.parse_args()

basedir = getenv('PANDA_FLATDIR') + ''
infile = basedir+'/%s.root'%(args.pd)
figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Tools.Load import Load
from PandaCore.Tools.Misc import *
from PandaCore.Tools.root_interface import Selector
from math import sqrt
Load('GraphAsymmErrDrawer')

import PandaAnalysis.SMH.BaselineSel as sel 

met = 'pfmet' if args.pd == 'SingleElectron' else 'pfUWmag'
base_cut = sel.cuts['singlemu'] if args.pd == 'MET' else sel.cuts['singleel']
met_trig = sel.triggers['met']

plot = root.GraphAsymmErrDrawer()
plot.SetTDRStyle()
plot.SetLumi(50)
plot.DrawEmpty(True)
#plot.SetAutoRange(False)
root.gStyle.SetOptStat(0)


s_mono = Selector()
s_di = Selector()
f = root.TFile.Open(infile); t = f.Get('events')
branches = list(set([met, 'pfmet', met_trig]))
s_mono.read_tree(t, branches = branches, cut = tAND(sel.monojet, base_cut))
s_di.read_tree(t, branches = branches, cut = tAND(sel.dijet, base_cut))

def eff(s, k):
  triggered = ( s[met_trig] != 0 )

  bins = [0, 50, 100, 150, 200, 250, 350, 500, 700, 1000]
  labels = {'pfmet':'E_{T}^{miss}', 'pfUWmag':'U(#mu)'}

  h_inc = s.draw(k, weight=None, vbins=bins)
  h_trg = s.draw(k, weight=None, vbins=bins, mask=(s[met_trig] != 0))

  h_ratio = root.TGraphAsymmErrors()
  h_ratio.BayesDivide(h_trg, h_inc)

  h_ratio.GetXaxis().SetTitle(labels[k])
  h_ratio.GetYaxis().SetTitle('Efficiency')
  h_ratio.SetMaximum(1.4)
  h_ratio.SetMinimum(0)

  return h_ratio


fout = root.TFile(getenv('CMSSW_BASE')+
                    '/src/PandaAnalysis/data/smh17/trig/met_trig_%s.root'%args.pd, 
                  'RECREATE')

for b in branches:
  if b == met_trig:
    continue 
  plot.Clear()
  plot.InitLegend(.6,.75)
  plot.AddCMSLabel()
  plot.AddLumiLabel()

  h_mono = eff(s_mono, b)
  h_di = eff(s_di, b)

  h_di.SetLineColor(root.kRed)

  plot.AddGraph(h_mono, 'p_{T}^{jet 1}>100', 1, 1, 'ep')
  plot.AddGraph(h_di, 'p_{T}^{jet 1}>60,p_{T}^{jet 2}>35', 2, 1, 'ep')

  plot.Draw(args.outdir+'/','eff_%s_%s'%(b, args.pd))

  fout.WriteTObject(h_mono, 'monojet_%s_%s'%(b, args.pd))
  fout.WriteTObject(h_di, 'dijet_%s_%s'%(b, args.pd))

fout.Close()
