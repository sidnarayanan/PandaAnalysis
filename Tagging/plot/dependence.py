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
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
from PandaCore.Tools.root_interface import Selector
Load('HistogramDrawer')

lumi=35800
plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.AddCMSLabel()
plot.AddSqrtSLabel()
#plot.SetAutoRange(False)
plot.InitLegend()
root.gStyle.SetOptStat(0)

nBins=5

f = root.TFile(args.infile)
t = f.Get('events')

h_effs = {}
h_errs = {}

base_base_cut = 'partonPt>300 && partonPt<470 && fabs(partonEta)<2.4'

if args.isfake:
  base_cut = tAND(base_base_cut, 'partonIsReco==1')
else:
  base_cut = tAND(base_base_cut, 'partonIsReco==1 && partonSize<1.44')


s = Selector()

s.read_tree(t, branches = ['partonPt','npv','1'], cut = base_base_cut)
h_inc = {'partonPt' : s.draw('partonPt', '1', fbins = (nBins, 300, 470)),
         'npv' : s.draw('npv', '1', fbins=(nBins,1,41))}

s.read_tree(t, branches = ['partonPt','npv','1'], cut = tAND(base_cut, 'top_ecf_bdt>0.43 && fj1MSD>110 && fj1MSD<210'))
h_bdt = {'partonPt' : s.draw('partonPt', '1', fbins = (nBins, 300, 470)),
         'npv' : s.draw('npv', '1', fbins=(nBins,1,41))}

s.read_tree(t, branches = ['partonPt','npv','1'], cut = tAND(base_cut, 'fj1Tau32SD<0.57 && fj1MSD>110 && fj1MSD<210'))
h_tau = {'partonPt' : s.draw('partonPt', '1', fbins = (nBins, 300, 470)),
         'npv' : s.draw('npv', '1', fbins=(nBins,1,41))}

if not args.isfake:
    s.read_tree(t, branches = ['partonPt','npv','1'], cut = base_cut)
    h_reco = {'partonPt' : s.draw('partonPt', '1', fbins = (nBins, 300, 470)),
             'npv' : s.draw('npv', '1', fbins=(nBins,1,41))}

for ib in xrange(1, nBins+1):
    for v in h_inc.values():
        v.SetBinError(ib,0)

for v, label in [('partonPt','gen p_{T} [GeV]'), ('npv','NPV')]:
    h_bdt[v].Divide(h_inc[v])
    h_tau[v].Divide(h_inc[v])
    if not args.isfake:
        h_reco[v].Divide(h_inc[v])
    plot.Reset(False)
    h_tau[v].GetXaxis().SetTitle(label)
    if not args.isfake:
        h_tau[v].SetMaximum(1)
    else:
        h_tau[v].SetMaximum(0.015)
    h_tau[v].GetYaxis().SetTitle('Efficiency')
    plot.AddHistogram(h_tau[v], '30% #tau_{32}^{SD}+m_{SD} WP', 11, root.kBlue, "e hist")
    plot.AddHistogram(h_bdt[v], '30% BDT+m_{SD} WP', 12, root.kRed, "e hist")
    if not args.isfake:
        plot.AddHistogram(h_reco[v], 'CA15 jet reco\'d', 10, root.kBlack, "e hist")
    if args.isfake:
        filename = 'qcd_'
    else:
        filename = 'sig_'
    filename += v 
    plot.Draw(args.outdir, filename)

