#!/usr/bin/env python

from os import system,getenv
from sys import argv
from PandaCore.Tools.script import * 
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
from PandaCore.Drawers.plot_utility import *

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
dataDir = baseDir

args = parse('--outdir', 
             ('--cut', {'default':'1==1'}),
             '--region')
lumi = 35800.
blind=True
region = args.region
sname = argv[0]

### DEFINE REGIONS ###

### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.Stack(False)
plot.SetTDRStyle()
plot.InitLegend(0.25, 0.65, 0.9, 0.88, 2)
plot.DrawMCErrors(True)
#plot.AddCMSLabel()
plot.cut = args.cut 
plot.SetEvtNum("eventNumber")
plot.SetLumi(lumi/1000)
plot.SetNormFactor(True)
plot.AddSqrtSLabel()
plot.do_overflow = False
plot.do_underflow = False
plot.n_threads = 10

weight = 'normalizedWeight * sf_qcdV * sf_ewkV' 
plot.mc_weight = weight

#logger.info('cut',plot.cut)
#logger.info('weight',plot.mc_weight)

#plot.add_systematic('QCD scale','scaleUp','scaleDown',root.kRed+2)
#plot.add_systematic('PDF','pdfUp','pdfDown',root.kBlue+2)

### DEFINE PROCESSES ###
znunu         = Process('Z(#nu#nu)+jets',root.kZjets)
zll           = Process('Z(ll)+jets',root.kExtra2)
wjets         = Process('W+jets',root.kWjets)
ttbar         = Process('t#bar{t}',root.kTTbar)
signal1        = Process('m_{V}=1.75 TeV, m_{#chi}=1 GeV',root.kSignal)
signal2        = Process('m_{#phi}=1.7 TeV, m_{#psi}=100 GeV',root.kSignal2)
processes = [ttbar,wjets, znunu, zll, signal1, signal2]

### ASSIGN FILES TO PROCESSES ###
zll.add_file(baseDir+'ZJets.root')
znunu.add_file(baseDir+'ZtoNuNu.root')
ttbar.add_file(baseDir+'TTbar.root')
wjets.add_file(baseDir+'WJets.root')
signal1.add_file(baseDir+'Vector_MonoTop_NLO_Mphi-1750_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph.root')
signal2.add_file(baseDir+'Scalar_MonoTop_LO_Mphi-1700_Mchi-100_13TeV-madgraph.root')

for p in processes:
    plot.add_process(p)


recoil=FDistribution("pfmet",200,1000,20,"p_{T}^{miss} [GeV]","Arbitrary units")
plot.add_distribution(recoil)

# plot.add_distribution(FDistribution('nJet',0.5,6.5,6,'N_{jet}','Events'))
# plot.add_distribution(FDistribution('npv',0,45,45,'N_{PV}','Events'))
plot.add_distribution(FDistribution('fjMSD[0]',0,250,20,'fatjet m_{SD} [GeV]','Arbitrary units'))
plot.add_distribution(FDistribution('fjPt[0]',250,1000,20,'fatjet p_{T} [GeV]','Arbitrary units'))
plot.add_distribution(FDistribution('top_ecf_bdt',-1,1,20,'Top BDT','Arbitrary units'))
# plot.add_distribution(FDistribution('fj1MaxCSV',0,1,20,'fatjet max CSV','Events'))
# plot.add_distribution(FDistribution('fj1Tau32',0,1,20,'fatjet #tau_{32}','Events'))
# plot.add_distribution(FDistribution('fj1Tau32SD',0,1,20,'fatjet #tau_{32}^{SD}','Events'))
# plot.add_distribution(FDistribution('jet1CSV',0,1,20,'jet 1 CSV','Events',filename='jet1CSV'))
# plot.add_distribution(FDistribution('dphipfmet',0,3.14,20,'min#Delta#phi(jet,E_{T}^{miss})','Events'))
#plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
plot.draw_all(args.outdir+'/'+region+'_')
