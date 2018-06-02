#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse

basedir = getenv('PANDA_FLATDIR')

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--infile',metavar='infile',type=str,default=basedir+'/TTbar.root')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--isfake',action='store_true')
args = parser.parse_args()

figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Utils.load import *
from PandaCore.Tools.Misc import *
from PandaCore.Tools.root_interface import Selector
Load('HistogramDrawer')

lumi=35800
plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.SetLumi(lumi/1000.)
plot.AddCMSLabel()
plot.AddLumiLabel()
#plot.SetAutoRange(False)
plot.InitLegend()
plot.AddPlotLabel('110 < m_{SD} < 210 GeV',.18,.77,False,42,.04)
root.gStyle.SetOptStat(0)

def get_max(h):
  m = 0
  for ib in xrange(1,h.GetNbinsX()+1): m = max(m,h.GetBinContent(ib))
  return m

nBins=10

f = root.TFile(args.infile)
t = f.Get('events')

h_effs = {}
h_errs = {}

if args.isfake:
  base_cut = 'fj1MSD>110 && fj1MSD<210'
else:
  base_cut = 'fj1MSD>110 && fj1MSD<210 && fj1IsMatched==1 && fj1GenSize<1.44'

s = Selector()

s.read_tree(t, branches = ['fj1Pt','normalizedWeight'], cut = base_cut)
h_inc = s.draw('fj1Pt', 'normalizedWeight', fbins = (nBins, 250, 1000))

colors = [root.kRed, root.kRed+3, root.kBlack]
counter=10
for c,col in zip([0.45,0.1,-0.5], colors):
  s.read_tree(t, branches = ['fj1Pt','normalizedWeight'], cut = tAND(base_cut, 'top_ecf_bdt>%f'%c))
  h_effs[c] = s.draw('fj1Pt', 'normalizedWeight', fbins = (nBins, 250, 1000))
  h_effs[c].Divide(h_inc)
  h_effs[c].SetMaximum(1.4*get_max(h_effs[c]))
  h_effs[c].GetXaxis().SetTitle('fatjet p_{T} [GeV]')
  h_effs[c].GetYaxis().SetTitle('Efficiency')
  h_effs[c].SetLineColor(col)
  plot.AddHistogram(h_effs[c], 'BDT > %.2f'%(c), counter, col, "hist")
  counter += 1

filename = 'pt_dep'
if args.isfake:
  filename = 'fake_'+filename
plot.Draw(args.outdir,filename)
