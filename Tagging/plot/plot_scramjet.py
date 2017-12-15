#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse
from collections import namedtuple

### SET GLOBAL VARIABLES ###
basedir = getenv('SCRAMJETFLAT')+'/' 
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--basedir',metavar='basedir',type=str,default=None)
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
args = parser.parse_args()

ylabel = 'Jets'

sname = argv[0]
if args.basedir:
    basedir = args.basedir

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
from PandaCore.Drawers.plot_utility import *
import PandaCore.Drawers.plot_utility as pu
import PandaAnalysis.Tagging.Selection as sel

pu.tree_name = 'puppiCA15'

### DEFINE REGIONS ###

### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.SetTDRStyle()
#plot.InitLegend()
plot.InitLegend(0.59,0.7,0.94,0.92)
plot.DrawMCErrors(True)
plot.AddCMSLabel(0.19,0.84,' Simulation')
plot.cut = 'mSD>110 && mSD<210 && pt<1000'
plot.SetNormFactor(True)
plot.AddSqrtSLabel()
plot.AddPlotLabel("110 < m_{SD} < 210 GeV",.19,.77,False,42,.04)

plot.mc_weight = 'normalizedWeight'

### DEFINE PROCESSES ###
matched = Process('Top quark jets',root.kExtra1)
matched.additional_cut = 'matched==1 && gensize<1.2'
matched.additional_weight = 'ptweight'

qcd = Process('q/g jets',root.kExtra2,1)
#qcd = Process('q/g jets',root.kExtra2,1)
qcd.dashed = True
qcd.additional_weight = 'ptweight_analytic'
processes = [qcd,matched]

### ASSIGN FILES TO PROCESSES ###
matched.add_file(basedir+'/ZpTT.root')
qcd.add_file(basedir+'/QCD.root')

for p in processes:
    plot.add_process(p)

plot.add_distribution(FDistribution('top_ecfv8_bdt',-1,1,50,'BDT output','Normalized distribution'))
plot.add_distribution(FDistribution('tau32SD',0,1,50,'#tau_{32}^{SD}','Normalized distribution'))


### DRAW AND CATALOGUE ###
plot.draw_all(args.outdir)
