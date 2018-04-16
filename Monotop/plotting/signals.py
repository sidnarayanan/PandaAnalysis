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

cut = 'pfmet>250 && nFatjet==1 && fj1Pt>250'
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

weight = 'normalizedWeight' 
plot.mc_weight = weight
plot.AddPlotLabel('m_{#chi}=1 GeV, g^{V}_{q}=0.25, g^{V}_{#chi}=1',.18,.77,False,42,.04)

#mVs = [300,500,1000,1500,2500]
mVs = [300,1000,2500]

### DEFINE PROCESSES ###
procs = [] 
counter=root.kExtra1
for m in mVs:
    matched = Process('m_{V}=%.1f TeV [top]'%(m/1000.),counter)
    matched.additional_cut = 'fj1IsMatched==1 && fj1GenSize<1.44'

    wmatched = Process('m_{V}=%.1f, [W]'%(m/1000.),counter)
    wmatched.additional_cut = '(fj1IsMatched==0 || fj1GenSize>1.44) && fj1IsWMatched==1 && fj1GenWSize<1.44'
    wmatched.dashed = True

    unmatched = Process('m_{V}=%.1f, [none]'%(m/1000.),counter)
    unmatched.additional_cut = '(fj1IsMatched==0 || fj1GenSize>1.44) && (fj1IsWMatched==0 || fj1GenWSize>1.44)'
    unmatched.dotted = True

    counter += 1
    for p in [unmatched,wmatched,matched]:
    #for p in [matched]:
        p.add_file(baseDir + '/Vector_MonoTop_NLO_Mphi-%i_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph.root'%m)
        procs.append(p)

for p in procs:
    plot.add_process(p)

plot.add_distribution(FDistribution('max(genAntiTopPt,genTopPt)',0,1000,40,'gen t p_{T} [GeV]','a.u.',filename='genTopPt'))
plot.add_distribution(FDistribution('fj1Pt',250,1000,40,'fatjet p_{T} [GeV]','a.u.'))
plot.add_distribution(FDistribution('fj1MSD',0,400,40,'fatjet m_{SD} [GeV]','a.u.'))
plot.add_distribution(FDistribution('top_ecf_bdt',-1,1,40,'Top BDT','a.u.'))

### DRAW AND CATALOGUE ###
plot.draw_all(args.outdir+'/')
