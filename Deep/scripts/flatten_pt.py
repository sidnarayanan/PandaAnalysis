#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse

basedir = getenv('PANDA_FLATDIR')

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--proc')
parser.add_argument('--plot',type=str,default=None)
args = parser.parse_args()

argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
from PandaCore.Tools.root_interface import Selector
if args.plot:
    Load('HistogramDrawer')

    plot = root.HistogramDrawer()
    plot.SetTDRStyle()
    plot.AddCMSLabel()
    plot.AddSqrtSLabel()
    plot.InitLegend()
    root.gStyle.SetOptStat(0)

binning = (40, 450, 1200)

n_partons = {
        'Top' : 3,
        'QCD' : 1,
        'Higgs' : 2,
        'W' : 2,
        }

f = root.TFile(basedir + '/' + args.proc + '.root')
t = f.Get('events')

s = Selector()

s.read_tree(t, branches = ['fj1RawPt'], cut = 'fj1RawPt>450 && fj1RawPt<1200')

h = s.draw('fj1RawPt', fbins = binning)

h_inv = h.Clone()
for ib in xrange(1, h_inv.GetNbinsX()+1):
    if h.GetBinContent(ib):
      h_inv.SetBinContent(ib, 1)
      h_inv.SetBinError(ib, 0)

h_inv.Divide(h)

for ib in xrange(1, h_inv.GetNbinsX()+1):
    # print ib, h.GetBinContent(ib), h_inv.GetBinContent(ib)
    if h.GetBinContent(ib) == 0:
      h_inv.SetBinContent(ib, 0)

fout = root.TFile.Open(getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/deep/flatten.root', 'UPDATE')
fout.WriteTObject(h_inv, 'h_'+args.proc, 'overwrite')
fout.Close()

h_inv.Scale(h.Integral()/100.)

fout = root.TFile.Open(getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/deep/flatten_scaled.root', 'UPDATE')
fout.WriteTObject(h_inv, 'h_'+args.proc, 'overwrite')
fout.Close()


if args.plot:
    plot.AddHistogram(h_inv, '', 11, root.kBlue, "e hist")
    plot.Draw(args.plot, 'weight_'+args.proc)
