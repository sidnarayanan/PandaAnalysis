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
plot.SetNormFactor(True)
#plot.SetAutoRange(False)
plot.InitLegend()
root.gStyle.SetOptStat(0)

nBins=15

f = root.TFile(args.infile)
t = f.Get('events')

base_cut = 'partonPt>300 && partonPt<470 && fabs(partonEta)<2.4 && partonIsReco==1'
if not args.isfake:
    base_cut = tAND(base_cut,'partonSize<1.44')


s = Selector()

s.read_tree(t, branches = ['fj1MSD','1'], cut = base_cut)
h_inc = s.draw('fj1MSD', '1', fbins = (nBins, 40, 340))

s.read_tree(t, branches = ['fj1MSD','1'], cut = tAND(base_cut, 'top_ecf_bdt>0.43'))
h_bdt = s.draw('fj1MSD', '1', fbins = (nBins, 40, 340))

s.read_tree(t, branches = ['fj1MSD','1'], cut = tAND(base_cut, 'fj1Tau32SD<0.43'))
h_tau = s.draw('fj1MSD', '1', fbins = (nBins, 40, 340))


v = 'fj1MSD'
label = 'm_{SD} [GeV]'

h_inc.GetXaxis().SetTitle(label)
h_inc.GetYaxis().SetTitle('Arbitrary units')
plot.AddHistogram(h_inc, 'Inclusive', 10, root.kBlack, "e hist")
plot.AddHistogram(h_tau, '30% #tau_{32}^{SD}+m_{SD} WP', 11, root.kBlue, "e hist")
plot.AddHistogram(h_bdt, '30% BDT+m_{SD} WP', 12, root.kRed, "e hist")
if args.isfake:
    filename = 'sculpt_qcd_'
else:
    filename = 'sculpt_sig_'
filename += v 
plot.Draw(args.outdir, filename)

