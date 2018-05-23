#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
dataDir = baseDir

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cat',metavar='cat',type=str,default='resolved')
args = parser.parse_args()

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.SMH.VBFSel as sel
from PandaCore.Drawers.plot_utility import *

cut = sel.cuts['signal']
boosted = 'fjPt>400 && fjDoubleCSV>0.9 && jot12Mass>2000'
if args.cat == 'resolved':
    cut = tAND(cut, tNOT(boosted))
else:
    cut = tAND(cut, boosted)

### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.Stack(False)
plot.SetTDRStyle()
plot.InitLegend()
plot.AddCMSLabel()
#plot.SetNormFactor(True)
plot.AddSqrtSLabel()
plot.cut = cut
#plot.do_overflow = True
plot.do_underflow = True

weight = sel.weights['signal']%1
plot.mc_weight = weight

### DEFINE PROCESSES ###
vbf           = Process('qq#rightarrowH#rightarrowbb',root.kSignal2)
qcd           = Process('QCD', root.kQCD)
tt           = Process('t#bar{t}', root.kTTbar)

### ASSIGN FILES TO PROCESSES ###
vbf.add_file(baseDir+'VBFHbb.root')
tt.add_file(baseDir+'TTbar.root')
for f in ['QCD_ht1000to1500.root','QCD_ht100to200.root','QCD_ht1500to2000.root',
          'QCD_ht2000toinf.root','QCD_ht200to300.root','QCD_ht300to500.root',
          'QCD_ht500to700.root','QCD_ht700to1000.root']:
    qcd.add_file(baseDir+f)

processes = [vbf,qcd,tt]

for p in processes:
    plot.add_process(p)

plot.add_distribution(FDistribution('jot12Mass',0,4000,20,'m_{jj} [GeV]','a.u.'))
plot.add_distribution(FDistribution('jot12DEta',0,10,20,'#Delta#eta(j_{1},j_{2})','a.u.'))
plot.add_distribution(FDistribution("fabs(jot12DPhi)",0,3.142,20,"#Delta #phi leading jets","a.u.",filename='jot12DPhi'))
plot.add_distribution(FDistribution("jotEta[0]",-5,5,20,"Jet 1 #eta","a.u."))
plot.add_distribution(FDistribution("jotEta[1]",-5,5,20,"Jet 2 #eta","a.u."))
plot.add_distribution(FDistribution("jotPt[0]",30,500,20,"Jet 1 p_{T} [GeV]","a.u."))
plot.add_distribution(FDistribution("jotPt[1]",30,500,20,"Jet 2 p_{T} [GeV]","a.u."))
plot.add_distribution(FDistribution("hbbm_dreg",0,300,30,"Resolved m_{H} [GeV]","a.u."))
if args.cat == 'boosted':
    plot.add_distribution(FDistribution("fjMSD_corr",0,250,30,"Boosted corrected m_{H} [GeV]","a.u."))
    plot.add_distribution(FDistribution("fjMSD",0,250,30,"Boosted m_{H} [GeV]","a.u."))

### DRAW AND CATALOGUE ###
plot.draw_all(args.outdir+'/'+args.cat+'_')
