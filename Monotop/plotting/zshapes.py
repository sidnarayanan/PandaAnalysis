#!/usr/bin/env python

from PandaCore.Tools.script import * 
import numpy as np

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
args = parse('--outdir', ('--cut', {'default':'1==1'}))
from PandaCore.Drawers.plot_utility import * 

### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.Logy(True)
plot.Stack(True)
#plot.DrawEmpty(True)
plot.Ratio(True)
plot.FixRatio(0.6)
plot.SetRatioLabel('NLO/LO')
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.do_overflow = False
plot.do_underflow = False

plot.mc_weight = 'normalizedWeight'
plot.cut = 'genJet1Pt>100 && fabs(genJet1Eta)<2.4 && genBosonPt>100 && lheHT>160'

plot.add_systematic('LO QCD uncertainty', 'scaleUp', 'scaleDown', 2)

### DEFINE PROCESSES ###
nlo     = Process('NLO (QCD) pp#rightarrowZ',root.kData); 
nlo.additional_weight = 'normalizedWeight' 
lo     = Process('LO pp#rightarrowZ',root.kZjets); 
processes = [lo, nlo]

### ASSIGN FILES TO PROCESSES ###
nlo.add_file(baseDir+'ZtoNuNu_nlo.root')
lo.add_file(baseDir+'ZtoNuNu.root')

for p in processes:
  plot.add_process(p)

ptbins = [160,200,250,300,350,400,450,500,550,600,700,800,900,1000,1100,1200,1400]
ptbins = np.linspace(160,1600,20)

plot.add_distribution(VDistribution('trueGenBosonPt',ptbins,'p_{T}^{Z} [GeV]','d#sigma/dp_{T} [pb/GeV]',filename='zpt'))
plot.add_distribution(FDistribution('genBosonMass',60,120,20,'m_{Z} [GeV]','d#sigma/dm_{Z} [pb/GeV]',filename='zm'))
plot.add_distribution(VDistribution('lheHT',ptbins,'H_{T} [GeV]', 'd#sigma/dH_{T} [pb/GeV]', filename='zht'))
plot.draw_all(args.outdir+'/')
