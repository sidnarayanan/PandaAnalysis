#!/usr/bin/env python

from os import system,getenv
from sys import argv
from PandaCore.Tools.script import * 
import PandaCore.Tools.Functions
import PandaAnalysis.Monotop.CombinedBVetoSelection as sel
from PandaCore.Drawers.plot_utility import *

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
dataDir = baseDir
args = parse('--outdir', '--region',
             ('--cut', {'default':'1==1'}),
             '--cat')
lumi = 36000.
blind=True
region = args.region


### DEFINE REGIONS ###
cut = tAND(sel.cuts[args.region],args.cut)
if args.cat:
    if args.cat == 'tight':
        cut = tAND('top_ecf_bdt>0.45', cut)
    else:
        cut = tAND('top_ecf_bdt<0.45&&top_ecf_bdt>0.1', cut)

### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.Stack(True)
plot.Ratio(True)
plot.FixRatio(0.5)
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
#plot.AddCMSLabel(0.15,0.94,'Preliminary')
#plot.AddCMSLabel(0.15,0.94,'')
plot.cut = cut
plot.SetEvtNum("eventNumber")
plot.SetLumi(lumi/1000)
plot.AddLumiLabel(True)
plot.do_overflow = True
plot.do_underflow = True
plot.n_threads = 5

weight = sel.weights[region]%lumi
plot.mc_weight = weight
if args.cat=='tight':
    plot.AddPlotLabel('BDT>0.45',.18,.8,False,42,.06)
elif args.cat=='loose':
    plot.AddPlotLabel('0.1<BDT<0.45',.18,.8,False,42,.06)

#plot.add_systematic('QCD scale','scaleUp','scaleDown',root.kRed+2)
#plot.add_systematic('PDF','pdfUp','pdfDown',root.kBlue+2)

### DEFINE PROCESSES ###
zjets         = Process('Z+jets',root.kZjets)
wjets         = Process('W+jets',root.kWjets)
diboson       = Process('Diboson',root.kDiboson)
ttbar         = Process('t#bar{t}',root.kTTbar)
ttg           = Process('t#bar{t}#gamma',root.kTTbar)
singletop     = Process('Single t',root.kST)
singletopg    = Process('t#gamma',root.kST)
qcd           = Process("QCD multijet",root.kQCD)
gjets         = Process('#gamma+jets',root.kGjets)
data          = Process("Data",root.kData)
signal        = Process('m_{V}=1.75 TeV, m_{#chi}=1 GeV',root.kSignal)
processes = [qcd,diboson,singletop,wjets,ttbar,zjets]
if 'dimuon' in region:
    processes = [diboson,singletop,ttbar,zjets]

### ASSIGN FILES TO PROCESSES ###
if 'signal' in region or 'qcd' in region:
    zjets.add_file(baseDir+'ZtoNuNu.root')
    signal.add_file(baseDir+'Vector_MonoTop_NLO_Mphi-1750_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph.root')
else:
    zjets.add_file(baseDir+'ZJets.root')
    #zjets.add_file(baseDir+'ZJets_nlo.root')
wjets.add_file(baseDir+'WJets.root')
diboson.add_file(baseDir+'Diboson.root')
ttbar.add_file(baseDir+'TTbar.root');
singletop.add_file(baseDir+'SingleTop.root')
ttg.add_file(baseDir+'TTbar_GJets.root');
singletopg.add_file(baseDir+'SingleTop_tG.root')
if 'pho' in region:
    processes = [qcd,singletopg,ttg,gjets]
    #processes = [qcd,gjets]
    gjets.add_file(baseDir+'GJets.root')
    qcd.add_file(baseDir+'SinglePhoton.root')
    qcd.additional_cut = sel.triggers['pho']
    qcd.use_common_weight = False
    qcd.additional_weight = 'sf_phoPurity'
else:
    qcd.add_file(baseDir+'QCD.root')

if any([x in region for x in ['singlemuonw','singleelectronw']]):
    processes = [qcd,diboson,singletop,zjets,ttbar,wjets,]
if any([x in region for x in ['singlemuontop','singleelectrontop']]):
    processes = [qcd,diboson,singletop,zjets,wjets,ttbar]
if any([x in region for x in ['signal','muon','qcd']]):
    data.additional_cut = sel.triggers['met']
    data.add_file(dataDir+'MET.root')
    lep='#mu'
    lepname = 'muon'
elif 'electron' in region:
    if 'di' in region:
        data.additional_cut = tOR(sel.triggers['ele'],sel.triggers['pho'])
    else:
        data.additional_cut = sel.triggers['ele']
    data.add_file(dataDir+'SingleElectron.root')
    lep='e'
    lepname = 'electron'
elif region=='photon':
    data.additional_cut = sel.triggers['pho']
    data.add_file(dataDir+'SinglePhoton.root')


processes.append(data)

for p in processes:
    plot.add_process(p)

recoilBins = [250,280,310,350,400,450,600,1000]
nRecoilBins = len(recoilBins)-1

### CHOOSE DISTRIBUTIONS, LABELS ###
if 'signal' in region or 'qcd' in region:
    recoil=VDistribution("pfmet",recoilBins,"p_{T}^{miss} [GeV]","Events/GeV")
elif any([x in region for x in ['singlemuonw','singleelectronw','singlemuontop','singleelectrontop','singlemuon','singleelectron']]):
    recoil=VDistribution("pfUWmag",recoilBins,"U(%s) [GeV]"%(lep),"Events/GeV")
    plot.add_distribution(FDistribution('%sPt[0]'%lepname,0,1000,20,'Leading %s p_{T} [GeV]'%lep,'Events/40 GeV'))
    plot.add_distribution(FDistribution('%sEta[0]'%lepname,-2.5,2.5,20,'Leading %s #eta'%lep,'Events/bin'))
elif any([x in region for x in ['dielectron','dimuon']]):
    recoil=VDistribution("pfUZmag",recoilBins,"U(%s%s) [GeV]"%(lep,lep),"Events/GeV")
    plot.add_distribution(FDistribution('diLepMass',60,120,20,'m_{ll} [GeV]','Events/3 GeV'))
    plot.add_distribution(FDistribution('%sPt[0]'%lepname,0,1000,20,'Leading %s p_{T} [GeV]'%lep,'Events/40 GeV'))
    plot.add_distribution(FDistribution('%sEta[0]'%lepname,-2.5,2.5,20,'Leading %s #eta'%lep,'Events/bin'))
    plot.add_distribution(FDistribution('%sPt[1]'%lepname,0,1000,20,'Subleading %s p_{T} [GeV]'%lep,'Events/40 GeV'))
    plot.add_distribution(FDistribution('%sEta[1]'%lepname,-2.5,2.5,20,'Subleading %s #eta'%lep,'Events/bin'))
elif region=='photon':
    recoil=VDistribution("pfUAmag",recoilBins,"U(#gamma) [GeV]","Events/GeV")
    plot.add_distribution(FDistribution('loosePho1Pt',0,1000,20,'Leading #gamma p_{T} [GeV]','Events/40 GeV'))
    plot.add_distribution(FDistribution('loosePho1Eta',-2.5,2.5,20,'Leading #gamma #eta','Events/bin'))

#recoil.calc_chi2 = True
plot.add_distribution(recoil)

plot.add_distribution(FDistribution('nJet',0.5,8.5,8,'N_{jet}','Events'))
plot.add_distribution(FDistribution('npv',0,45,45,'N_{PV}','Events'))
plot.add_distribution(FDistribution('fjMSD[0]',50,250,10,'fatjet m_{SD} [GeV]','Events'))
plot.add_distribution(FDistribution('fjPt[0]',200,1000,20,'fatjet p_{T} [GeV]','Events'))
plot.add_distribution(FDistribution('top_ecf_bdt',-1,1,20,'BDT Output','Events/0.1 Units'))
plot.add_distribution(FDistribution('fjMaxCSV[0]',0,1,20,'fatjet max CSV','Events'))
plot.add_distribution(FDistribution('fjTau[0]32',0,1,20,'fatjet #tau_{32}','Events'))
plot.add_distribution(FDistribution('fjTau[0]32SD',0,1,20,'fatjet #tau_{32}^{SD}','Events'))
#plot.add_distribution(FDistribution('jet1CSV',0,1,20,'jet 1 CSV','Events',filename='jet1CSV'))
plot.add_distribution(FDistribution('dphipfmet',0,3.14,20,'min#Delta#phi(jet,E_{T}^{miss})','Events'))
plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
if args.cat:
    region += '_' + args.cat
plot.draw_all(args.outdir+'/'+region+'_')
