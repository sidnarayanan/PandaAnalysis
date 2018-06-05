#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--wp',type=str,default='tight')
parser.add_argument('--pt',type=int,default=-1)
args = parser.parse_args()
lumi = 36000.
sname = argv[0]
wp = 0.45 if args.wp == 'tight' else 0.1 

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.Tagging.TnPSel as sel 
from PandaCore.Drawers.plot_utility import *

### DEFINE REGIONS ###
cut = tAND(sel.cuts['tag'], sel.triggers['mu'])
if args.pt >= 0:
    cut = tAND(cut, sel.pt_cut(args.pt))
weight = sel.weights['tag']

### LOAD PLOTTING UTILITY ###

plot = PlotUtility()
plot.Stack(True)
plot.Ratio(True)
plot.FixRatio(1)
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.SetLumi(lumi/1000)
plot.AddLumiLabel(True)
plot.do_overflow = False
plot.do_underflow = False
plot.mc_weight = weight%lumi
plot.n_threads = 1

plot.cut = cut

### DEFINE PROCESSES ###
p1 = Process('1-prong', root.kWjets)
p2 = Process('2-prong', root.kExtra1)
p3 = Process('3-prong', root.kTTbar)
p2.additional_cut = '(fjIsMatched==0||fjGenSize>1.44)'
p3.additional_cut = '(fjIsMatched==1&&fjGenSize<1.44)'
data = Process('Data', root.kData)

baseDir += '/tnp_1.44/'
data.add_file(baseDir+'Data.root')
p1.add_file(baseDir+'1.root')
p2.add_file(baseDir+'2.root')
p3.add_file(baseDir+'3.root')
#data.add_file(baseDir+'/SingleMuon.root')
#p1.add_file(baseDir+'/WJets.root')
#p2.add_file(baseDir+'/Diboson.root')
#p2.add_file(baseDir+'/SingleTop.root')
#p2.add_file(baseDir+'/TTbar.root')
#p3.add_file(baseDir+'/SingleTop.root')
#p3.add_file(baseDir+'/TTbar.root')

for p in [p1,p2,p3,data]:
    plot.add_process(p)

plot.add_distribution(FDistribution('fjMSD',50,350,30,'CA15 m_{SD} [GeV]','Events/10 GeV'))

fname = '/CAT_%s_'%args.wp
if args.pt >= 0:
    fname += '%i_'%args.pt 

# first pass
plot.cut = tAND(cut, 'top_ecf_bdt>%f'%wp)
plot.draw_all(args.outdir+fname.replace('CAT','pass'))

# then fail
plot.cut = tAND(cut, 'top_ecf_bdt<%f'%wp)
plot.draw_all(args.outdir+fname.replace('CAT','fail'))
