#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse

basedir = getenv('PANDA_FLATDIR')
submit_name = getenv('SUBMIT_NAME')

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--proc')
parser.add_argument('--plot',type=str,default=None)
args = parser.parse_args()

argv=[]

import ROOT as root
from PandaCore.Utils.load import *
from PandaCore.Tools.Misc import *
from PandaCore.Tools.root_interface import Selector
if args.plot:
    Load('HistogramDrawer')

    plot = root.HistogramDrawer()
    plot.SetTDRStyle()
    plot.AddSqrtSLabel()
    plot.InitLegend()
    root.gStyle.SetOptStat(0)

binning = (40, 450, 1200)

n_partons = {
        'Top_lo' : 3,
        'Top' : 3,
        'QCD' : 1,
        'Higgs' : 2,
        'W' : 2,
        }
try:
    n_parton = n_partons[args.proc]
except KeyError:
    n_parton = n_partons[args.proc.split('_')[0]]

f = root.TFile(basedir + '/' + args.proc + '.root')
t = f.Get('events')

s = Selector()

s.read_tree(t, branches = ['genFatJetPt'], cut = 'genFatJetPt>450 && genFatJetPt<1200 && genFatJetNProngs==%i'%n_parton)

h = s.draw('genFatJetPt', fbins = binning)

h_inv = h.Clone()
for ib in xrange(1, h_inv.GetNbinsX()+1):
    if h.GetBinContent(ib):
      h_inv.SetBinContent(ib, 1)
      h_inv.SetBinError(ib, 0)

h_inv.Divide(h)
h_inv.GetXaxis().SetTitle('jet p_{T}')
h_inv.GetYaxis().SetTitle('jet weight')

for ib in xrange(1, h_inv.GetNbinsX()+1):
    # print ib, h.GetBinContent(ib), h_inv.GetBinContent(ib)
    if h.GetBinContent(ib) == 0:
        h_inv.SetBinContent(ib, 0)

print 'Saving histogram as h_%s_%s'%(args.proc, submit_name)

fout = root.TFile.Open(getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/deep/flatten_gen.root', 'UPDATE')
fout.WriteTObject(h_inv, 'h_%s_%s'%(args.proc, submit_name), 'overwrite')
fout.Close()

h_inv.Scale(h.Integral()/50.)

fout = root.TFile.Open(getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/deep/flatten_gen_scaled.root', 'UPDATE')
fout.WriteTObject(h_inv, 'h_%s_%s'%(args.proc, submit_name), 'overwrite')
fout.Close()


if args.plot:
    plot.AddHistogram(h_inv, '', 11, root.kBlue, "e hist")
    plot.Draw(args.plot, 'weight_'+args.proc)
