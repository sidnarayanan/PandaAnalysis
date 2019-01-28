#!/usr/bin/env python

from os import system,getenv
from PandaCore.Tools.script import * 

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
dataDir = baseDir

args = parse('--outdir', ('--cat', {'default':'loose'}))
lumi = 36000.
blind=True
region = 'signal' 

import PandaCore.Tools.Functions
import PandaAnalysis.VBF.PandaSelection as sel
from PandaCore.Drawers.plot_utility import *

### DEFINE REGIONS ###

cut = sel.cuts[region]
cut = cut.replace('jotPt[0]>80 && jotPt[1]>40', 'jotPt[0]>0 && jotPt[1]>0')
if args.cat!='loose':
    cut = tAND(cut,getattr(sel,args.cat))
cut = 'pfmet>250 && dphipfmet>0.5 && jotEta[0]*jotEta[1]<0 '
#cut = 'jotPt[0]>80 && jotPt[1]>40'

### LOAD PLOTTING UTILITY ###
BLIND=True

plot = PlotUtility()
plot.Stack(False)
plot.SetLumi(lumi/1000)
plot.DrawMCErrors(True)
plot.SetTDRStyle()
#plot.InitLegend(0.2,0.75,0.9,0.9,3)
# plot.AddCMSLabel()
plot.SetNormFactor(True)
plot.cut = cut
print cut 
plot.AddLumiLabel(True)
plot.do_overflow = False 
plot.do_underflow = False
plot.SetAbsMin(0.0001)

#weight = sel.weights[region].replace('sf_l1','1')%lumi
#weight = 'normalizedWeight'
weight = '1'
plot.mc_weight = weight

### DEFINE PROCESSES ###
zjets         = Process('Z+jets [QCD]',root.kZjets) # ; zjets.dashed = True
wjets         = Process('W+jets [QCD]',root.kWjets) # ; wjets.dashed = True
zjets_ewk     = Process('Z+jets [EW]',root.kExtra4) # ; zjets_ewk.dotted = True
wjets_ewk     = Process('W+jets [EW]',root.kExtra5) # ; wjets_ewk.dotted = True
vbf           = Process('qq#rightarrowqqH(inv)',root.kSignal)
ggf           = Process('ggF H(inv)',root.kSignal2)

### ASSIGN FILES TO PROCESSES ###
zjets.add_file(baseDir+'ZtoNuNu.root')
zjets_ewk.add_file(baseDir+'ZtoNuNu_EWK.root')
wjets.add_file(baseDir+'WJets.root')
wjets_ewk.add_file(baseDir+'WJets_EWK.root')
vbf.add_file(baseDir+'vbfHinv_mAll.root')
ggf.add_file(baseDir+'ggFHinv_mAll.root')

processes = [ggf, vbf]
#processes = [zjets, wjets, zjets_ewk, wjets_ewk, vbf]

for p in processes:
    plot.add_process(p)

recoilBins = np.array([200., 230., 260.0, 290.0, 320.0, 350.0, 390.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0])
recoilBins /= 1000
nRecoilBins = len(recoilBins)-1

recoil=VDistribution("pfmet/1000",recoilBins,"PF MET [TeV]","Events/TeV",filename='pfmet')
#plot.add_distribution(recoil)

plot.add_distribution(FDistribution('jot12Mass/1000',0,5,15,'m_{jj} [TeV]','Arbitrary units',filename='jot12Mass'))
plot.add_distribution(FDistribution('jot12DEta',0,8,15,'#Delta#eta_{jj}','Arbitrary units'))
plot.add_distribution(FDistribution("fabs(jot12DPhi)",0,3.142,15,"#Delta#phi_{jj}","Arbitrary units",filename='jot12DPhi'))
plot.add_distribution(FDistribution("fabs(jotEta[0])",0,5,15,"Jet 1 |#eta|","Arbitrary units"))
plot.add_distribution(FDistribution("fabs(jotEta[1])",0,5,15,"Jet 2 |#eta|","Arbitrary units"))
plot.add_distribution(FDistribution("jotPt[0]",0,700,15,"Jet 1 p_{T} [GeV]","Arbitrary units"))
plot.add_distribution(FDistribution("jotPt[1]",0,500,15,"Jet 2 p_{T} [GeV]","Arbitrary units"))
plot.add_distribution(FDistribution("jotE[0]",0,2000,15,"Jet 1 E [GeV]","Arbitrary units"))
plot.add_distribution(FDistribution("jotE[1]",0,2000,15,"Jet 2 E [GeV]","Arbitrary units"))
#plot.add_distribution(FDistribution("nJot-0.5",1.5,6.5,5,"N_{jet}","Arbitrary units",filename='nJot'))
#plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
region = '%s/%s'%(args.cat,region)
plot.draw_all(args.outdir+'/'+region+'_')
