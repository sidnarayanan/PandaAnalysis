#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
dataDir = baseDir

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--bin',metavar='bin',type=int,default=None)
args = parser.parse_args()
lumi = 35900.
sname = argv[0]

mjj_binning = [200, 400, 600, 900, 1200, 1500, 2000, 2750, 3500, 5000]
mjj_lo = mjj_binning[args.bin]
mjj_hi = mjj_binning[args.bin+1]

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.VBF.PandaSelection as sel
from PandaCore.Drawers.plot_utility import *

### DEFINE REGIONS ###

cut = tAND(sel.cuts['qcd'], sel.mjj)
cut = tAND(cut, 'jot12Mass>%i && jot12Mass<%i'%(mjj_lo, mjj_hi))

### LOAD PLOTTING UTILITY ###

plot = PlotUtility()
plot.Stack(True)
plot.SetTDRStyle()
# plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.AddPlotLabel('%i < m_{jj} < %i GeV'%(mjj_lo, mjj_hi),
                  0.18, 0.8, False, 42, 0.05)
plot.cut = cut
plot.AddSqrtSLabel()
plot.do_overflow = True
plot.do_underflow = True

plot.mc_weight = sel.weights['qcd']%lumi

### DEFINE PROCESSES ###
qcd           = Process("QCD",root.kQCD)

### ASSIGN FILES TO PROCESSES ###
qcd.add_file(baseDir+'QCD.root')
plot.add_process(qcd)

# recoilBins = range(0, 200, 20) + [200., 230., 260.0, 290.0, 320.0, 350.0, 390.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]
recoilBins = [0, 50, 100, 150, 200, 300, 400, 600, 800, 1200]
nRecoilBins = len(recoilBins)-1

### CHOOSE DISTRIBUTIONS, LABELS ###
recoil=VDistribution("pfmet",recoilBins,"PF MET [GeV]","Events/GeV")
plot.add_distribution(recoil)

### DRAW AND CATALOGUE ###
plot.draw_all(args.outdir+'/mjj_bin%i'%args.bin+'_')
