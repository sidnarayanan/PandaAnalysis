#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange
import numpy as np

basedir = getenv('PANDA_FLATDIR') + ''


parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--nmu',metavar='nmu',type=int,default=1)
args = parser.parse_args()

if args.nmu == 1:
  infile = basedir+'/WJets.root'
  datafile = basedir + '/SingleMuon.root'
else:
  infile = basedir+'/ZJets.root'
  datafile = basedir + '/SingleMuon.root'
figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Tools.Load import Load
from PandaCore.Tools.Misc import *
from PandaCore.Tools.root_interface import Selector
from math import sqrt
Load('GraphAsymmErrDrawer')
root.gROOT.SetBatch()

lumi=36000
import PandaAnalysis.VBF.PandaSelection as sel 
if args.nmu == 1:
  base_cut = sel.cuts['singlemuon']
  metnomu = 'pfUWmag'
else:
  base_cut = sel.cuts['dimuon']
  metnomu = 'pfUZmag'
base_cut = tAND(base_cut, '(trigger&8)!=0')

plot = root.GraphAsymmErrDrawer()
plot.SetTDRStyle()
plot.SetLumi(lumi/1000.)
plot.DrawEmpty(True)
#plot.SetAutoRange(False)
root.gStyle.SetOptStat(0)


s = Selector()
f = root.TFile.Open(infile); t = f.Get('events')
branches = ['trigger&1',  'barrelHTMiss', 'trigger&8','normalizedWeight']
s.read_tree(t, branches = branches, cut = base_cut)

sdata = Selector()
fdata = root.TFile.Open(datafile); tdata = fdata.Get('events')
branches.remove('normalizedWeight')
sdata.read_tree(tdata, branches = branches, cut = base_cut)

triggered = ( s['trigger&1'] != 0 )
triggered_data = ( sdata['trigger&1'] != 0 )

ratios_to_plot = {
  'barrelHTMiss' : ([0, 80, 120, 160, 200, 240, 280, 320, 360, 400, 450, 500, 600, 1000],'Barrel H_{T}^{miss} [GeV]'),
}

for k,v in ratios_to_plot.iteritems():
    plot.Clear()
    plot.InitLegend()
    plot.AddCMSLabel()
    plot.AddLumiLabel()

    h_inc = s.draw(k, weight='normalizedWeight', vbins=v[0])
    h_trig = s.draw(k, weight='normalizedWeight', mask=triggered, vbins=v[0])

    hdata_inc = sdata.draw(k, weight=None, vbins=v[0])
    hdata_trig = sdata.draw(k, weight=None, mask=triggered_data, vbins=v[0])

    h_ratio_inc = root.TGraphAsymmErrors()
    h_ratio_inc.BayesDivide(h_trig, h_inc)
    h_ratiodata = root.TGraphAsymmErrors()
    h_ratiodata.BayesDivide(hdata_trig, hdata_inc)
    h_ratio_inc.GetXaxis().SetTitle(v[1])
    h_ratio_inc.GetYaxis().SetTitle('E_{T}^{miss} trig eff')
    h_ratio_inc.SetMaximum(1.4)
    h_ratio_inc.SetMinimum(0)

    plot.AddGraph(h_ratio_inc, 'MC', 1, 1, 'ep')
    plot.AddGraph(h_ratiodata, 'Data', 2, 1, 'ep')

    plot.Draw(args.outdir,'compare_trigeff_nmu%i'%args.nmu+k)
