#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange
import numpy as np

basedir = getenv('PANDA_FLATDIR') + '/' 

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--proc',metavar='proc',type=str)
parser.add_argument('--outdir',metavar='outdir',type=str)
args = parser.parse_args()

figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Tools.Load import Load
from PandaCore.Tools.Misc import *
from PandaCore.Tools.root_interface import Selector
from math import sqrt
Load('HistogramDrawer')
root.gROOT.SetBatch()

lumi=36000
import PandaAnalysis.VBF.PandaSelection as sel 

plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.Logy()
plot.DrawEmpty(True)

#plotr = root.HistogramDrawer()
#plot.SetRatioStyle()
#plot.DrawEmpty(True)

s = Selector()
s_nlo = Selector()

cut = 'trueGenBosonPt>160'
if args.proc == 'WJets':
    cut = tAND(cut, 'genBosonMass<70 && genBosonMass<90')
else:
    cut = tAND(cut, 'genBosonMass<80 && genBosonMass<100')

infile = basedir + args.proc + '.root'
f = root.TFile.Open(infile); t = f.Get('events')
branches = ['normalizedWeight', 'trueGenBosonPt', 'genMjj', 'genJet1Pt', 'genJet2Pt']
s.read_tree(t, branches = branches, cut = cut)

infile = basedir + args.proc + '_nlo.root'
f_nlo = root.TFile.Open(infile); t_nlo = f_nlo.Get('events')
s_nlo.read_tree(t_nlo, branches = branches, cut = cut)

def draw(x, xbins,  xlabel, ylabel):
    hlo = s.draw(x, weight='normalizedWeight', vbins=xbins)
    hnlo = s_nlo.draw(x, weight='normalizedWeight', vbins=xbins)
    hlo.Scale(1,'width'); hnlo.Scale(1,'width')
    
    plot.Reset()
    plot.InitLegend()
    plot.AddCMSLabel()
    plot.AddSqrtSLabel()
    hnlo.GetXaxis().SetTitle(xlabel)
    hnlo.GetYaxis().SetTitle(ylabel)
    plot.AddHistogram(hnlo, 'NLO', 7)
    plot.AddHistogram(hlo, 'LO', 8)
    plot.Draw(args.outdir,x+'_dist_'+args.proc)

#    hratio = hnlo.Clone(); hratio.Divide(hlo)
#    plotr.Reset()
#    plotr.AddCMSLabel()
#    plotr.AddSqrtSLabel()
#    hratio.GetXaxis().SetTitle(xlabel)
#    hratio.GetYaxis().SetTitle('NLO/LO')
#    plotr.AddHistogram(hratio, 'NLO/LO', 7, -1, 'hist e')
#    plotr.Draw(args.outdir,x+'_ratio_'+args.proc)
  
 

toplot = [
    ('trueGenBosonPt', [160,180,200,250,300,400,500,600,700,800,1000,1200,1400,2000], 'p_{T}^{V} [GeV]', 'd#sigma/dp_{T}^{V}'), 
    ('genJet1Pt', np.arange(0,1000,40), 'p_{T}^{jet 1} [GeV]', 'd#sigma/dp_{T}^{jet 1}'), 
    ('genJet2Pt', np.arange(0,1000,40), 'p_{T}^{jet 2} [GeV]', 'd#sigma/dp_{T}^{jet 2}'), 
    ('genMjj', [0,200,400,800,1000,1400,1800,2400,3000,4000,5000], 'm_{jj} [GeV]', 'd#sigma/dm_{jj}'), 
    ]

for x in toplot:
    draw(*x)
