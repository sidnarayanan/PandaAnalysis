#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
dataDir = baseDir

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
parser.add_argument('--region',metavar='region',type=str,default='signal')
parser.add_argument('--cat',type=str,default='mjj')
parser.add_argument('--cr',type=str,default='lowdphi'),
args = parser.parse_args()
lumi = 36000.
blind=True
region = args.region
sname = argv[0]
do_jec = True

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.VBF.PandaSelection as sel
#import PandaAnalysis.VBF.TriggerSelection as sel
from PandaCore.Drawers.plot_utility import *

### DEFINE REGIONS ###

cut = tAND(sel.cuts[args.region],args.cut)
cut = cut.replace('pfmet>250','pfmet<250 && pfmet>100')
if args.cat!='loose':
    cut = tAND(cut,getattr(sel,args.cat))
if args.cr == 'lowdphi':
    cut = cut.replace('dphipfmet>0.5','dphipfmet<0.1')

### LOAD PLOTTING UTILITY ###
BLIND=True

plot = PlotUtility()
plot.Stack(True)
plot.Ratio(True)
plot.FixRatio(1)
plot.SetAbsMin(0.001)
plot.SetTDRStyle("vbf")
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.cut = cut
plot.SetEvtNum("eventNumber")
plot.SetLumi(lumi/1000)
plot.AddLumiLabel(True)
plot.do_overflow = True
plot.do_underflow = True

if args.cat == 'cnc':
    weight = sel.weights_cnc[region]%lumi
else:
    weight = sel.weights[region]%lumi
plot.mc_weight = tTIMES('0.5',weight)

### DEFINE PROCESSES ###
zjets         = Process('Z+jets [QCD]',root.kZjets)
wjets         = Process('W+jets [QCD]',root.kWjets)
zjets_ewk     = Process('Z+jets [EWK]',root.kExtra3)
wjets_ewk     = Process('W+jets [EWK]',root.kExtra2)
qcd           = Process("QCD",root.kQCD)
data          = Process("Data",root.kData)

### ASSIGN FILES TO PROCESSES ###
zjets.add_file(baseDir+'ZtoNuNu.root')
zjets_ewk.add_file(baseDir+'ZtoNuNu_EWK.root')
data.add_file(baseDir+'JetHT.root')
#data.additional_cut = sel.triggers['met']
wjets.add_file(baseDir+'WJets.root')
wjets_ewk.add_file(baseDir+'WJets_EWK.root')
qcd.add_file(baseDir+'QCD.root')

processes = [wjets_ewk,wjets,zjets_ewk,zjets,qcd,data]

for p in processes:
    plot.add_process(p)


recoilBins = [200., 230., 260.0, 290.0, 320.0, 350.0, 390.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]
nRecoilBins = len(recoilBins)-1


plot.add_distribution(VDistribution('jot12Mass',[200, 600, 1200, 2000, 2750, 5000],'m_{jj} [GeV]','Events/GeV'))
plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
region = args.cr
region = '%s/%s'%(args.cat,region)
plot.draw_all(args.outdir+'/'+region+'_')
