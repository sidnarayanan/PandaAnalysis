#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange

basedir = getenv('PANDA_FLATDIR')
infile = basedir+'/QCD.root'


parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
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
cut = 'jot1Eta*jot2Eta<0 && jot1Pt>80 && jot2Pt>40 && fabs(jot1Eta)<4.7 && fabs(jot2Eta)<4.7 && nTau==0 && jetNMBtags==0 && jot1VBFID==1 && dphipfmet>0.5'
weight = 'normalizedWeight'



s = Selector()
f = root.TFile.Open(infile); t = f.Get('events')
branches = ['pfmet', 'jot12Mass', 
            weight]
s.read_tree(t, branches = branches, cut = cut)

met_bins = array('f',[50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 220, 240, 260, 280, 300, 340, 380, 420, 460, 500, 550, 600, 650, 700, 750, 800, 900, 1000])
mjj_bins = array('f',[200, 400, 600, 900, 1200, 1500, 2000, 2750, 3500, 5000])
h = root.TH2F('h','h',len(mjj_bins)-1, mjj_bins, len(met_bins)-1, met_bins)
h.GetXaxis().SetTitle('m_{jj} [GeV]')
h.GetYaxis().SetTitle('E_{T}^{miss} [GeV]')
h = s.draw(fields=['jot12Mass', 'pfmet'], weight=weight, hbase=h)

plot = root.CanvasDrawer()
plot.SetTDRStyle()
root.gStyle.SetPadRightMargin(0.15)
c = root.TCanvas()
c.cd() 
c.SetLogx(); c.SetLogy(); c.SetLogz()
plot.AddCMSLabel(.16,.94)
plot.AddSqrtSLabel()
plot.SetCanvas(c)
root.gStyle.SetOptStat(0)
root.gStyle.SetNumberContours(999);
root.gStyle.SetPalette(root.kBird)

h.Draw('colz')

plot.Draw(args.outdir,'qcd_met_mjj')

fout = root.TFile.Open(args.outdir+'/qcd_met_mjj.root', 'RECREATE')
fout.WriteTObject(h, 'h_mjj_met')
fout.Close()

