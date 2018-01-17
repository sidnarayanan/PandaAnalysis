#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
dataDir = baseDir

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
args = parser.parse_args()
lumi = 36000.
blind=True
region = 'signal'
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

scales = {

        'Diboson' : "(0.893678*(200.0<jot12Mass&&jot12Mass<400.0))+(0.893678*(400.0<jot12Mass&&jot12Mass<600.0))+(0.893678*(600.0<jot12Mass&&jot12Mass<900.0))+(0.893678*(900.0<jot12Mass&&jot12Mass<1200.0))+(0.893678*(1200.0<jot12Mass&&jot12Mass<1500.0))+(0.893678*(1500.0<jot12Mass&&jot12Mass<2000.0))+(0.893678*(2000.0<jot12Mass&&jot12Mass<2750.0))+(0.893678*(2750.0<jot12Mass&&jot12Mass<3500.0))+(0.893678*(3500.0<jot12Mass&&jot12Mass<5000.0))",
        'Top' : "(0.856676*(200.0<jot12Mass&&jot12Mass<400.0))+(0.856677*(400.0<jot12Mass&&jot12Mass<600.0))+(0.856676*(600.0<jot12Mass&&jot12Mass<900.0))+(0.856676*(900.0<jot12Mass&&jot12Mass<1200.0))+(0.856676*(1200.0<jot12Mass&&jot12Mass<1500.0))+(0.856676*(1500.0<jot12Mass&&jot12Mass<2000.0))+(0.856676*(2000.0<jot12Mass&&jot12Mass<2750.0))+(0.856676*(2750.0<jot12Mass&&jot12Mass<3500.0))+(0.856676*(3500.0<jot12Mass&&jot12Mass<5000.0))",
        'Z+jets [EWK]' : "(1.007564*(200.0<jot12Mass&&jot12Mass<400.0))+(1.025401*(400.0<jot12Mass&&jot12Mass<600.0))+(1.054702*(600.0<jot12Mass&&jot12Mass<900.0))+(1.076168*(900.0<jot12Mass&&jot12Mass<1200.0))+(1.102413*(1200.0<jot12Mass&&jot12Mass<1500.0))+(1.022041*(1500.0<jot12Mass&&jot12Mass<2000.0))+(0.977125*(2000.0<jot12Mass&&jot12Mass<2750.0))+(0.827909*(2750.0<jot12Mass&&jot12Mass<3500.0))+(0.733391*(3500.0<jot12Mass&&jot12Mass<5000.0))",
        'W+jets [EWK]' : "(1.058023*(200.0<jot12Mass&&jot12Mass<400.0))+(1.085203*(400.0<jot12Mass&&jot12Mass<600.0))+(1.107709*(600.0<jot12Mass&&jot12Mass<900.0))+(1.145824*(900.0<jot12Mass&&jot12Mass<1200.0))+(1.152217*(1200.0<jot12Mass&&jot12Mass<1500.0))+(1.094762*(1500.0<jot12Mass&&jot12Mass<2000.0))+(1.023695*(2000.0<jot12Mass&&jot12Mass<2750.0))+(0.875813*(2750.0<jot12Mass&&jot12Mass<3500.0))+(0.827997*(3500.0<jot12Mass&&jot12Mass<5000.0))",
        'QCD' : "(0.954892*(200.0<jot12Mass&&jot12Mass<400.0))+(0.928491*(400.0<jot12Mass&&jot12Mass<600.0))+(0.940957*(600.0<jot12Mass&&jot12Mass<900.0))+(0.930812*(900.0<jot12Mass&&jot12Mass<1200.0))+(0.932788*(1200.0<jot12Mass&&jot12Mass<1500.0))+(0.931599*(1500.0<jot12Mass&&jot12Mass<2000.0))+(0.936846*(2000.0<jot12Mass&&jot12Mass<2750.0))+(0.934052*(2750.0<jot12Mass&&jot12Mass<3500.0))+(0.914625*(3500.0<jot12Mass&&jot12Mass<5000.0))",
        'Z+jets [QCD]' : "(1.007564*(200.0<jot12Mass&&jot12Mass<400.0))+(1.025401*(400.0<jot12Mass&&jot12Mass<600.0))+(1.054702*(600.0<jot12Mass&&jot12Mass<900.0))+(1.076168*(900.0<jot12Mass&&jot12Mass<1200.0))+(1.102413*(1200.0<jot12Mass&&jot12Mass<1500.0))+(1.022041*(1500.0<jot12Mass&&jot12Mass<2000.0))+(0.977125*(2000.0<jot12Mass&&jot12Mass<2750.0))+(0.827909*(2750.0<jot12Mass&&jot12Mass<3500.0))+(0.733391*(3500.0<jot12Mass&&jot12Mass<5000.0))",
        'W+jets [QCD]' : "(1.008800*(200.0<jot12Mass&&jot12Mass<400.0))+(1.045349*(400.0<jot12Mass&&jot12Mass<600.0))+(1.053005*(600.0<jot12Mass&&jot12Mass<900.0))+(1.098035*(900.0<jot12Mass&&jot12Mass<1200.0))+(1.099920*(1200.0<jot12Mass&&jot12Mass<1500.0))+(1.044082*(1500.0<jot12Mass&&jot12Mass<2000.0))+(0.979727*(2000.0<jot12Mass&&jot12Mass<2750.0))+(0.830397*(2750.0<jot12Mass&&jot12Mass<3500.0))+(0.765759*(3500.0<jot12Mass&&jot12Mass<5000.0))",


        }

### DEFINE REGIONS ###

cut = tAND(tAND(sel.cuts[region],sel.mjj),'200<jot12Mass && jot12Mass<5000')

### LOAD PLOTTING UTILITY ###
BLIND=False

plot = PlotUtility()
plot.Stack(True)
plot.Ratio(True)
plot.FixRatio(0.5)
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

weight = sel.weights[region]%lumi
plot.mc_weight = weight

### DEFINE PROCESSES ###
zjets         = Process('Z+jets [QCD]',root.kZjets)
wjets         = Process('W+jets [QCD]',root.kWjets)
zjets_ewk     = Process('Z+jets [EWK]',root.kExtra3)
wjets_ewk     = Process('W+jets [EWK]',root.kExtra2)
diboson       = Process('Diboson',root.kDiboson)
top           = Process('Top',root.kTTbar)
qcd           = Process("QCD",root.kQCD)
data          = Process("Data",root.kData)
#vbf           = Process('VBF H(inv)',root.kExtra4)

### ASSIGN FILES TO PROCESSES ###
zjets.add_file(baseDir+'ZtoNuNu.root')
zjets_ewk.add_file(baseDir+'ZtoNuNu_EWK.root')
data.add_file(baseDir+'MET.root')
data.additional_cut = sel.triggers['met']
wjets.add_file(baseDir+'WJets.root')
wjets_ewk.add_file(baseDir+'WJets_EWK.root')
diboson.add_file(baseDir+'Diboson.root')
top.add_file(baseDir+'TTbar.root');
top.add_file(baseDir+'SingleTop.root');
qcd.add_file(baseDir+'QCD.root')

processes = [qcd,diboson,top,wjets_ewk,zjets_ewk,wjets,zjets]
for p in processes:
    p.additional_weight = scales[p.name]
processes.append(data)

for p in processes:
    plot.add_process(p)

recoilBins = [200., 230., 260.0, 290.0, 320.0, 350.0, 390.0, 430.0, 470.0, 510.0, 550.0, 590.0, 640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0]
nRecoilBins = len(recoilBins)-1

if 'electron' in region:
    lep = 'e'
else:
    lep = '#mu'

### CHOOSE DISTRIBUTIONS, LABELS ###
'''
'''
if 'signal' in region:
    recoil=VDistribution("pfmet",recoilBins,"PF MET [GeV]","Events/GeV")
    plot.add_distribution(FDistribution('dphipfmet',0.5,3.2,20,'min #Delta#phi(U,jets)','Events'))
elif 'single' in region:
    recoil=VDistribution("pfUWmag",recoilBins,"PF U(%s) [GeV]"%(lep),"Events/GeV")
    plot.add_distribution(FDistribution('looseLep1Pt',0,1000,20,'Leading %s p_{T} [GeV]'%lep,'Events/40 GeV'))
    plot.add_distribution(FDistribution('looseLep1Eta',-2.5,2.5,20,'Leading %s #eta'%lep,'Events/bin'))
#    plot.add_distribution(FDistribution('dphipfUW',0.5,3.2,20,'min #Delta#phi(U,jets)','Events'))
elif any([x in region for x in ['dielectron','dimuon']]):
    recoil=VDistribution("pfUZmag",recoilBins,"PF U(%s%s) [GeV]"%(lep,lep),"Events/GeV")
    plot.add_distribution(FDistribution('diLepMass',60,120,20,'m_{ll} [GeV]','Events/3 GeV'))
    plot.add_distribution(FDistribution('looseLep1Pt',0,1000,20,'Leading %s p_{T} [GeV]'%lep,'Events/40 GeV'))
    plot.add_distribution(FDistribution('looseLep1Eta',-2.5,2.5,20,'Leading %s #eta'%lep,'Events/bin'))
    plot.add_distribution(FDistribution('looseLep2Pt',0,1000,20,'Subleading %s p_{T} [GeV]'%lep,'Events/40 GeV'))
    plot.add_distribution(FDistribution('looseLep2Eta',-2.5,2.5,20,'Subleading %s #eta'%lep,'Events/bin'))
#    plot.add_distribution(FDistribution('dphipfUZ',0.5,3.2,20,'min #Delta#phi(U,jets)','Events'))

plot.add_distribution(recoil)

plot.add_distribution(FDistribution('barrelHT',0,1000,20,'Barrel H_{T} [GeV]','Events/50 GeV'))
plot.add_distribution(FDistribution('barrelHTMiss',0,1000,20,'Barrel H_{T}^{miss} [GeV]','Events/50 GeV'))
#plot.add_distribution(FDistribution('barrelJet1Pt',0,1000,20,'Barrel jet 1 p_{T} [GeV]','Events/50 GeV'))
plot.add_distribution(FDistribution('jot12Mass',0,4000,20,'m_{jj} [GeV]','Events/200 GeV'))
plot.add_distribution(FDistribution('jot12DEta',0,4,20,'#Delta#eta(j_{1},j_{2})','Events'))
#plot.add_distribution(FDistribution('jot12DEta',0,8.5,20,'#Delta#eta(j_{1},j_{2})','Events'))
plot.add_distribution(FDistribution("fabs(jot12DPhi)",0,2,20,"#Delta #phi leading jets","Events",filename='jot12DPhi'))
plot.add_distribution(FDistribution("TMath::Sqrt(jot12DPhi*jot12DPhi + jot12DEta*jot12DEta)",0,8,20,"#Delta R leading jets","Events",filename='jot12DR'))
#plot.add_distribution(FDistribution("jot1CHF",0,1,20,"Jet 1 CHF","Events"))
#plot.add_distribution(FDistribution("jot1NHF",0,1,20,"Jet 1 NHF","Events"))
#plot.add_distribution(FDistribution("jot2CHF",0,1,20,"Jet 2 CHF","Events"))
#plot.add_distribution(FDistribution("jot2NHF",0,1,20,"Jet 2 NHF","Events"))
plot.add_distribution(FDistribution("jot1Eta",-5,5,20,"Jet 1 #eta","Events"))
plot.add_distribution(FDistribution("jot2Eta",-5,5,20,"Jet 2 #eta","Events"))
plot.add_distribution(FDistribution("jot1Phi",-5,5,20,"Jet 1 #phi","Events"))
plot.add_distribution(FDistribution("jot2Phi",-5,5,20,"Jet 2 #phi","Events"))
plot.add_distribution(FDistribution("jot1Pt",80,500,20,"Jet 1 p_{T} [GeV]","Events"))
plot.add_distribution(FDistribution("jot2Pt",40,500,20,"Jet 2 p_{T} [GeV]","Events"))
plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))
'''
'''

### DRAW AND CATALOGUE ###
plot.draw_all(args.outdir)
