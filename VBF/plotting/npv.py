#!/usr/bin/env python

from PandaCore.Tools.script import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.VBF.PandaSelection as sel
#import PandaAnalysis.VBF.TriggerSelection as sel
from PandaCore.Drawers.plot_utility import *

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
dataDir = baseDir

args = parse('--outdir', ('--cut', {'default':'1==1'}),
             '--region', 
             '--f16', '--f18', '--f17')
lumi = 36000.
region = args.region

### DEFINE REGIONS ###

cut = tAND(sel.cuts[args.region],args.cut)

### LOAD PLOTTING UTILITY ###
BLIND=False

plot = PlotUtility()
plot.SetTDRStyle()
plot.SetNormFactor(True)
#plot.InitLegend(0,0,0,0)
#plot.AddCMSLabel()
plot.cut = 'npv>0'
plot.SetLumi(lumi/1000)
plot.AddSqrtSLabel()
plot.do_overflow = True
plot.do_underflow = True

plot.mc_weight = '1'

### DEFINE PROCESSES ###
years = [16, 17, 18]
for i,y in enumerate(years):
     data = Process("20%i"%y,root.kData, 2*i)
     data.add_file(getattr(args, 'f%i'%y))
     plot.add_process(data)

plot.add_distribution(FDistribution('npv',0,60,60,'N_{PV}','Events'))

### DRAW AND CATALOGUE ###
plot.draw_all(args.outdir+'/comparison_')
