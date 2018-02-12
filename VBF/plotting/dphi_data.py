#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
dataDir = baseDir

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--region',metavar='region',type=str,default=None)
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
from PandaCore.Drawers.plot_utility import *

### DEFINE REGIONS ###

cut = tAND(sel.cuts[args.region],sel.mjj)
cut = tAND(cut, sel.triggers['met'])

### LOAD PLOTTING UTILITY ###

plot = PlotUtility()
plot.Ratio(True)
plot.FixRatio(0.5)
plot.SetTDRStyle("vbf")
plot.InitLegend()
plot.DrawMCErrors(False)
plot.AddCMSLabel()
plot.cut = cut
plot.SetLumi(lumi/1000)
plot.AddLumiLabel(True)
plot.do_overflow = True
plot.do_underflow = True
plot.SetNormFactor(True)
plot.SetRatioLabel('[0.2,0.5] / Sideband')

### DEFINE PROCESSES ###
data_in       = Process("0.2 < #Delta#phi < 0.5",root.kData)
data_in.additional_cut = 'fabs(jot12DPhi)>0.2 && fabs(jot12DPhi)<0.5'
data_lo       = Process("#Delta#phi < 0.2",root.kData,3)
data_lo.additional_cut = 'fabs(jot12DPhi)<0.2'
data_lo.ratio = True
data_out       = Process("0.5 < #Delta#phi < 0.8",root.kData,2)
data_out.additional_cut = 'fabs(jot12DPhi)>0.5 && fabs(jot12DPhi)<0.8'
data_out.ratio = True

### ASSIGN FILES TO PROCESSES ###
for d in [data_in, data_out, data_lo]:
    d.add_file(baseDir+'/MET.root')
    plot.add_process(d)

recoilBins = [200., 230., 260.0, 290.0, 320.0, 350.0, 390.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]
nRecoilBins = len(recoilBins)-1

### CHOOSE DISTRIBUTIONS, LABELS ###
'''
'''
recoil=VDistribution("pfmet",recoilBins,"PF MET [GeV]","a.u./GeV")
plot.add_distribution(recoil)
plot.add_distribution(FDistribution('pfmetphi',-3.2,3.2,20,'PF MET #phi','a.u.'))
plot.add_distribution(FDistribution('calomet/pfmet',0.4,1.2,20,'Calo MET / PF MET','a.u.',filename='caloOverPF'))
#
#
plot.add_distribution(FDistribution('barrelHT',0,1000,20,'Barrel H_{T} [GeV]','a.u./50 GeV'))
plot.add_distribution(FDistribution('barrelHTMiss',0,1000,20,'Barrel H_{T}^{miss} [GeV]','a.u./50 GeV'))
plot.add_distribution(FDistribution('barrelJet1Pt',0,1000,20,'Barrel jet 1 p_{T} [GeV]','a.u./50 GeV'))
plot.add_distribution(FDistribution('jot12Mass',0,4000,20,'m_{jj} [GeV]','a.u./200 GeV'))
plot.add_distribution(FDistribution('jot12DEta',0,4,20,'#Delta#eta(j_{1},j_{2})','a.u.'))
plot.add_distribution(FDistribution('dphipfmet',0,3.2,25,'min #Delta#phi(U,jets)','a.u.'))
plot.add_distribution(FDistribution("fabs(jot12DPhi)",0,3,30,"#Delta #phi leading jets","a.u.",filename='jot12DPhi'))
plot.add_distribution(FDistribution("TMath::Sqrt(jot12DPhi*jot12DPhi + jot12DEta*jot12DEta)",0,8,20,"#Delta R leading jets","a.u.",filename='jot12DR'))
#plot.add_distribution(FDistribution("jot1CHF",0,1,20,"Jet 1 CHF","a.u."))
#plot.add_distribution(FDistribution("jot1NHF",0,1,20,"Jet 1 NHF","a.u."))
#plot.add_distribution(FDistribution("jot2CHF",0,1,20,"Jet 2 CHF","a.u."))
#plot.add_distribution(FDistribution("jot2NHF",0,1,20,"Jet 2 NHF","a.u."))
plot.add_distribution(FDistribution("jot1Eta",-5,5,20,"Jet 1 #eta","a.u."))
plot.add_distribution(FDistribution("jot2Eta",-5,5,20,"Jet 2 #eta","a.u."))
plot.add_distribution(FDistribution("jot1Phi",-5,5,20,"Jet 1 #phi","a.u."))
plot.add_distribution(FDistribution("jot2Phi",-5,5,20,"Jet 2 #phi","a.u."))
plot.add_distribution(FDistribution("jot1Pt",80,500,20,"Jet 1 p_{T} [GeV]","a.u."))
plot.add_distribution(FDistribution("jot2Pt",40,500,20,"Jet 2 p_{T} [GeV]","a.u."))
plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
plot.draw_all(args.outdir+'/'+region+'_')
