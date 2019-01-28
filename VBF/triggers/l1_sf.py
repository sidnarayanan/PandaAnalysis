#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = '/data/t3home000/zdemirag/forSid/panda/v_004_16/'
dataDir = baseDir

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
parser.add_argument('--region',metavar='region',type=str,default='signal')
parser.add_argument('--cat',type=str,default='mjj')
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
import PandaAnalysis.VBF.TestSelection as sel
#import PandaAnalysis.VBF.TriggerSelection as sel
from PandaCore.Drawers.plot_utility import *

### DEFINE REGIONS ###

cut = tAND(sel.cuts[args.region],args.cut)
if args.cat!='loose':
    cut = tAND(cut,getattr(sel,args.cat))

### LOAD PLOTTING UTILITY ###

plot = PlotUtility()
plot.Stack(True)
plot.Ratio(True)
plot.FixRatio(0.4)
plot.SetTDRStyle("vbf")
plot.InitLegend(0.4)
plot.DrawMCErrors(False)
#plot.AddCMSLabel()
plot.cut = cut
plot.AddSqrtSLabel()
plot.do_overflow = True
plot.do_underflow = True

if args.cat == 'cnc':
    weight = sel.weights_cnc[region]%lumi
else:
    weight = sel.weights[region]%lumi
plot.mc_weight = '1' # weight

### DEFINE PROCESSES ###
vbf           = Process('Without trigger effects',root.kExtra1)
vbfl1           = Process('With trigger inefficiencies',root.kData)
vbfl1.dashed = True
vbfl1.additional_weight = 'sf_l1finor*sf_metTrigVBF'

### ASSIGN FILES TO PROCESSES ###
for p in [vbf,vbfl1]:
    p.add_file(baseDir+'vbfHinv_m125.root')
    plot.add_process(p)

recoilBins = [200., 230., 260.0, 290.0, 320.0, 350.0, 390.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]
nRecoilBins = len(recoilBins)-1

if 'electron' in region:
    lep = 'e'
else:
    lep = '#mu'

### CHOOSE DISTRIBUTIONS, LABELS ###
#recoil=VDistribution("pfmet",recoilBins,"PF MET [GeV]","Events/GeV")
#plot.add_distribution(FDistribution('dphipfmet',0.5,3.2,20,'min #Delta#phi(U,jets)','Events'))
#plot.add_distribution(recoil)

#plot.add_distribution(FDistribution('barrelHT',0,1000,20,'Barrel H_{T} [GeV]','Events/50 GeV'))
#plot.add_distribution(FDistribution('barrelJet1Pt',0,1000,20,'Barrel jet 1 p_{T} [GeV]','Events/50 GeV'))
plot.add_distribution(FDistribution('jot12Mass',200,4200,15,'m_{jj} [GeV]','Events/bin'))
#plot.add_distribution(FDistribution('jot12DEta',0,10,20,'#Delta#eta(j_{1},j_{2})','Events'))
#plot.add_distribution(FDistribution("fabs(jot12DPhi)",0,3.142,20,"#Delta #phi leading jets","Events",filename='jot12DPhi'))
#plot.add_distribution(FDistribution("jot1Eta",-5,5,20,"Jet 1 #eta","Events"))
#plot.add_distribution(FDistribution("jot2Eta",-5,5,20,"Jet 2 #eta","Events"))
#plot.add_distribution(FDistribution("jot1Pt",80,500,20,"Jet 1 p_{T} [GeV]","Events"))
#plot.add_distribution(FDistribution("jot2Pt",40,500,20,"Jet 2 p_{T} [GeV]","Events"))
#plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
region = '%s/%s'%(args.cat,region)
plot.draw_all(args.outdir+'/'+region+'_')
