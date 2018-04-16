#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--basedir',metavar='basedir',type=str,default=None)
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
args = parser.parse_args()
lumi = 36560.
sname = argv[0]
if args.basedir:
    baseDir = args.basedir

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
from PandaCore.Drawers.plot_utility import *

### DEFINE REGIONS ###

cut = '(pfmet>250||pfUWmag>250) && nFatjet==1 && fj1Pt>250'
cut = tAND(cut,args.cut)


### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.SetTDRStyle()
plot.InitLegend()
plot.AddCMSLabel()
plot.cut = cut
plot.SetEvtNum("eventNumber")
plot.SetLumi(lumi/1000)
plot.AddLumiLabel(True)
plot.do_overflow = True
plot.do_underflow = True
plot.SetNormFactor(True)

weight = '%f*normalizedWeight'%lumi
plot.mc_weight = weight


### DEFINE PROCESSES ###
procs = [] 
ttbar = Process('t#bar{t} [heavy]',root.kTTbar); ttbar.additional_cut = 'jet1Flav>0'
lfttbar = Process('t#bar{t} [light]',root.kTTbar); lfttbar.additional_cut = 'jet1Flav==0'; lfttbar.dashed=True
wjets = Process('W+jets [heavy]',root.kWjets); wjets.additional_cut = 'jet1Flav>0'
lfwjets = Process('W+jets [light]',root.kWjets); lfwjets.additional_cut = 'jet1Flav==0'; lfwjets.dashed=True

ttbar.add_file(baseDir+'TTbar.root')
lfttbar.add_file(baseDir+'TTbar.root')
wjets.add_file(baseDir+'WJets.root')
lfwjets.add_file(baseDir+'WJets.root')

procs = [ttbar,lfttbar,wjets,lfwjets]
procs.reverse()

for p in procs:
    plot.add_process(p)

plot.add_distribution(FDistribution('jet1Pt',0,1000,40,'Leading jet p_{T} [GeV]','a.u.'))
plot.add_distribution(FDistribution('fj1Pt',250,1000,40,'fatjet p_{T} [GeV]','a.u.'))

### DRAW AND CATALOGUE ###
plot.draw_all(args.outdir+'/')
