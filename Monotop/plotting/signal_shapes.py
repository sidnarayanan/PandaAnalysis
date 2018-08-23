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
lumi = 36000
sname = argv[0]
if args.basedir:
    baseDir = args.basedir

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
from PandaCore.Drawers.plot_utility import *

### DEFINE REGIONS ###

cut = 'pfmet>250 && nFatjet==1 && fj1Pt>250 && top_ecf_bdt > 0.1'
cut = tAND(cut,args.cut)


### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.SetTDRStyle()
plot.InitLegend()
#plot.AddCMSLabel(0.18,0.85,' Simulation')
plot.cut = cut
plot.AddSqrtSLabel()
plot.do_overflow = False
plot.do_underflow = False
plot.SetNormFactor(True)

weight = '%f*normalizedWeight' % lumi
plot.mc_weight = weight

SCALAR = False

if SCALAR:
    plot.AddPlotLabel('#splitline{m_{#psi}=100 GeV}{a_{q}=b_{q}=0.1, a_{#psi}=b_{#psi}=0.2}',.18,.81,False,42,.04)
    mVs = [2300,2900,3500,4100]
else:
    plot.AddPlotLabel('m_{#chi}=1 GeV, g^{V}_{q}=0.25, g^{V}_{#chi}=1',.18,.83,False,42,.04)
    mVs = [500,1000,1500,2500]

# mVs = [300,1000,2500]

### DEFINE PROCESSES ###
procs = [] 
counter=root.kExtra1
for m in mVs:
    if SCALAR:
        p = Process('m_{#phi}=%.1f TeV '%(m/1000.),counter)
        p.add_file(baseDir + '/Scalar_MonoTop_LO_Mphi-%i_Mchi-100_13TeV-madgraph.root'%m)
    else:
        p = Process('m_{V}=%.1f TeV '%(m/1000.),counter)
        p.add_file(baseDir + '/Vector_MonoTop_NLO_Mphi-%i_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph.root'%m)
    procs.append(p)
    counter += 1

for p in procs:
    plot.add_process(p)

if SCALAR:
    plot.add_distribution(FDistribution('pfmet',250,2500,25,'p_{T}^{miss} [GeV]','Arbitrary units'))
else:
    plot.add_distribution(FDistribution('pfmet',250,1000,25,'p_{T}^{miss} [GeV]','Arbitrary units'))

### DRAW AND CATALOGUE ###
if SCALAR:
    plot.draw_all(args.outdir+'/res_')
else:
    plot.draw_all(args.outdir+'/fcnc_')
