#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
dataDir = baseDir

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cat',type=str,default='loose')
args = parser.parse_args()
lumi = 35800.
blind=True
region = 'signal' 
sname = argv[0]

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.VBF.PandaSelection as sel
from PandaCore.Drawers.plot_utility import *

### DEFINE REGIONS ###

cut = sel.cuts[region]
for pt in ['jot1Pt','jot2Pt']:
    cut = removeCut(cut,pt)
cut = tAND('jot1Pt>0 && jot2Pt>0', cut)
if args.cat!='loose':
    cut = tAND(cut,getattr(sel,args.cat))

### LOAD PLOTTING UTILITY ###
BLIND=True

plot = PlotUtility()
plot.Stack(False)
plot.SetLumi(lumi/1000)
plot.SetTDRStyle()
plot.InitLegend()
plot.AddCMSLabel()
plot.SetNormFactor(True)
plot.cut = cut
plot.AddLumiLabel(True)
plot.do_overflow = True
plot.do_underflow = True

weight = sel.weights[region]%lumi
plot.mc_weight = weight

### DEFINE PROCESSES ###
zjets         = Process('Z+jets [QCD]',root.kZjets)
wjets         = Process('W+jets [QCD]',root.kWjets)
zjets_ewk     = Process('Z+jets [EWK]',root.kExtra3)
wjets_ewk     = Process('W+jets [EWK]',root.kExtra2)
vbf           = Process('VBF H(inv)',root.kSignal)
ggf           = Process('ggF H(inv)',root.kSignal2)

### ASSIGN FILES TO PROCESSES ###
zjets.add_file(baseDir+'ZtoNuNu.root')
zjets_ewk.add_file(baseDir+'ZtoNuNu_EWK.root')
wjets.add_file(baseDir+'WJets.root')
wjets_ewk.add_file(baseDir+'WJets_EWK.root')
vbf.add_file(baseDir+'vbfHinv_m125.root')
ggf.add_file(baseDir+'ggFHinv_m125.root')

processes = [zjets, wjets, zjets_ewk, wjets_ewk, ggf, vbf]

for p in processes:
    plot.add_process(p)

recoilBins = [200., 230., 260.0, 290.0, 320.0, 350.0, 390.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]
nRecoilBins = len(recoilBins)-1

recoil=VDistribution("pfmet",recoilBins,"PF MET [GeV]","a.u.")
plot.add_distribution(recoil)

plot.add_distribution(FDistribution('jot12Mass',0,4000,20,'m_{jj} [GeV]','a.u.'))
plot.add_distribution(FDistribution('jot12DEta',0,10,20,'#Delta#eta(j_{1},j_{2})','a.u.'))
plot.add_distribution(FDistribution("fabs(jot12DPhi)",0,3.142,20,"#Delta #phi leading jets","a.u.",filename='jot12DPhi'))
plot.add_distribution(FDistribution("jot1Eta",-5,5,20,"Jet 1 #eta","a.u."))
plot.add_distribution(FDistribution("jot2Eta",-5,5,20,"Jet 2 #eta","a.u."))
plot.add_distribution(FDistribution("jot1Pt",30,500,20,"Jet 1 p_{T} [GeV]","a.u."))
plot.add_distribution(FDistribution("jot2Pt",30,500,20,"Jet 2 p_{T} [GeV]","a.u."))
plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
region = '%s/%s'%(args.cat,region)
plot.draw_all(args.outdir+'/'+region+'_')
