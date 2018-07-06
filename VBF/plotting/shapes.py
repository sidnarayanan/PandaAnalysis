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
plot.InitLegend(0.64,0.55,0.92,0.9)
plot.AddCMSLabel(0.18,0.85,' Supplementary')
plot.SetNormFactor(True)
plot.cut = cut
plot.AddSqrtSLabel()
plot.do_overflow = True
plot.do_underflow = True

weight = sel.weights[region]%lumi
plot.mc_weight = weight

### DEFINE PROCESSES ###
vjets         = Process('V+jets [QCD]',root.kZjets, root.kRed)
vjets_ewk     = Process('V+jets [EW]',root.kZjets, root.kBlue)
vbf           = Process('qqH#rightarrowinv',root.kZjets, root.kBlack)
ggf           = Process('ggH#rightarrowinv',root.kZjets, root.kBlack); ggf.dashed = True

### ASSIGN FILES TO PROCESSES ###
vjets.add_file(baseDir+'ZtoNuNu.root')
vjets_ewk.add_file(baseDir+'ZtoNuNu_EWK.root')
vjets.add_file(baseDir+'WJets.root')
vjets_ewk.add_file(baseDir+'WJets_EWK.root')
vbf.add_file(baseDir+'vbfHinv_m125.root')
ggf.add_file(baseDir+'ggFHinv_m125.root')

processes = [vjets, vjets_ewk, vbf, ggf]

for p in processes:
    plot.add_process(p)

recoilBins = [200., 230., 260.0, 290.0, 320.0, 350.0, 390.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]
nRecoilBins = len(recoilBins)-1

recoil=VDistribution("pfmet",recoilBins,"PF MET [GeV]","Arbitrary units")
plot.add_distribution(recoil)

mjjbins = [0,200,400,600,800,1200,1600,2000,2500,3000,4000]
plot.add_distribution(VDistribution('jot12Mass',mjjbins,'m_{jj} [GeV]','Arbitrary units'))

plot.add_distribution(FDistribution('jot12DEta',0,8,16,'#Delta#eta(j_{1},j_{2})','Arbitrary units'))
plot.add_distribution(FDistribution("fabs(jot12DPhi)",0,3.142,12,"#Delta #phi leading jets","Arbitrary units",filename='jot12DPhi'))
plot.add_distribution(FDistribution("jot1Eta",-5,5,20,"Jet 1 #eta","Arbitrary units"))
plot.add_distribution(FDistribution("jot2Eta",-5,5,20,"Jet 2 #eta","Arbitrary units"))
plot.add_distribution(FDistribution("jot1Pt",30,500,20,"Jet 1 p_{T} [GeV]","Arbitrary units"))
plot.add_distribution(FDistribution("jot2Pt",30,500,20,"Jet 2 p_{T} [GeV]","Arbitrary units"))
#plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
region = '%s/%s'%(args.cat,region)
plot.draw_all(args.outdir+'/'+region+'_')
