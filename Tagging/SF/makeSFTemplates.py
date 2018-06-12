#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv,exit
import argparse

basedir = getenv('PANDA_FLATDIR')
figsdir = basedir+'/figs'

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--indir',metavar='indir',type=str,default=basedir)
parser.add_argument('--outdir',metavar='outdir',type=str,default=figsdir)
parser.add_argument('--disc',metavar='disc',type=str,default='top_ecfv8_bdt')
parser.add_argument('--tagged',metavar='tagged',type=str,default='True')
parser.add_argument('--sel',metavar='sel',type=str,default='tag')
parser.add_argument('--cut',metavar='cut',type=float,default=0.5)
parser.add_argument('--pt',metavar='pt',type=str,default='inc')
parser.add_argument('--shower',metavar='shower',type=str,default='pythia')
args = parser.parse_args()

figsdir = args.outdir+'/'+args.pt
basedir = args.indir
argv=[]

import ROOT as root
from PandaCore.Utils.load import *
from PandaCore.Tools.Misc import *
tcut = root.TCut
Load('Drawers','PlotUtility')

### SET GLOBAL VARIABLES ###
lumi = 12918.
logy=False
nlo = 'sf_ewkV*sf_qcdV'
if args.sel=='mistag':
  cut = 'nFatjet==1 && fj1Pt>250 && fj1MaxCSV<0.46 && nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && UWmag>250 && isojetNBtags==0'
  weight = '%f*normalizedWeight*sf_pu*sf_lep*%s*sf_sjbtag0*sf_btag0*sf_tt*sf_metTrig'%(lumi,nlo)
  label = 'mistag_'
elif args.sel=='photon':
  cut = 'nFatjet==1 && fj1Pt>250 && nLooseLep==0 && nLoosePhoton==1 && loosePho1IsTight==1 && nTau==0 && UAmag>250'
  weight = '%f*normalizedWeight*sf_pu*sf_lep*%s*sf_tt*sf_phoTrig*0.93'%(lumi,nlo)
  label = 'photon_'
else:
  cut = 'nFatjet==1 && fj1Pt>250 && fj1MaxCSV>0.46 && nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && UWmag>250 && isojetNBtags==1'
  weight = '%f*normalizedWeight*sf_pu*sf_lep*%s*sf_sjbtag1*sf_btag1*sf_tt*sf_metTrig'%(lumi,nlo)
  label = 'tag_'

disclabels = {
    'fj1Tau32SD':'#tau_{32}^{SD}',
    'top_ecfv6_bdt':'ECF BDT',
    'top_ecfv7_bdt':'ECF+#tau_{32}^{SD} BDT',
    'top_ecfv8_bdt':'ECF+#tau_{32}^{SD}+f_{rec} BDT',
    'top_ecf_bdt':'ECF+#tau_{32}^{SD}+f_{rec} BDT',
    }

label += args.disc+'_'
if args.tagged.upper()=='TRUE':
  label += 'pass_'
  if 'Tau' in args.disc:
    cut = tAND(cut,'%s<%f'%(args.disc,args.cut))
    plotlabel = '%s<%.2f'%(disclabels[args.disc],args.cut)
  else:
    cut = tAND(cut,'%s>%f'%(args.disc,args.cut))
    plotlabel = '%s>%.2f'%(disclabels[args.disc],args.cut)
else:
  label += 'fail_'
  if 'Tau' in args.disc:
    cut = tAND(cut,'%s>%f'%(args.disc,args.cut))
    plotlabel = '%s>%.2f'%(disclabels[args.disc],args.cut)
  else:
    cut = tAND(cut,'%s<%f'%(args.disc,args.cut))
    plotlabel = '%s<%.2f'%(disclabels[args.disc],args.cut)

if args.pt=='' or args.pt=='inc':
  cut = tAND(cut,'fj1Pt>250 && fj1Pt<1000')
  plotlabel = '#splitline{%s}{250 < p_{T} < 1000 GeV}'%plotlabel
elif args.pt=='lo':
  cut = tAND(cut,'fj1Pt>250 && fj1Pt<450')
  plotlabel = '#splitline{%s}{250 < p_{T} < 450 GeV}'%plotlabel
elif args.pt=='med':
  cut = tAND(cut,'fj1Pt>450 && fj1Pt<1000')
  plotlabel = '#splitline{%s}{450 < p_{T} < 500 GeV}'%plotlabel
elif args.pt=='hi':
  cut = tAND(cut,'fj1Pt>500 && fj1Pt<1000')
  plotlabel = '#splitline{%s}{500 < p_{T} < 1000 GeV}'%plotlabel

### LOAD PLOTTING UTILITY ###
plot = root.PlotUtility()
plot.Ratio(True) 
plot.SetTDRStyle()
plot.Stack(True)
plot.InitLegend()
plot.SetCut(tcut(cut))
plot.FixRatio(.5)
plot.SetLumi(lumi/1000)
plot.DrawMCErrors(True)
plot.SetNormFactor(False)
plot.AddCMSLabel()
plot.AddLumiLabel()
if plotlabel:
  plot.AddPlotLabel(plotlabel,.18,.77,False,42,.04)
plot.SetMCWeight(weight)

### DEFINE PROCESSES ###
prong1     = root.Process('1-prong',root.kWjets) 
prong2     = root.Process('2-prong',root.kExtra1); prong2.additionalCut = root.TCut("(fj1IsMatched==0||fj1GenSize>1.44)");
prong3     = root.Process('3-prong',root.kTTbar); prong3. additionalCut = root.TCut("(fj1IsMatched==1&&fj1GenSize<1.44)")
photon     = root.Process('#gamma+jets',root.kGjets)
qcd        = root.Process('QCD',root.kQCD)
data       = root.Process('Data',root.kData)
if args.sel=='photon':
  data.additionalCut = root.TCut('(trigger&4)!=0')
  processes = [photon,qcd,data]
else:
  data.additionalCut = root.TCut('(trigger&1)!=0')
  processes = [prong1,prong2,prong3,data]

### ASSIGN FILES TO PROCESSES ###
if args.sel=='photon':
  data.AddFile(basedir+'SinglePhoton.root')
  photon.AddFile(basedir+'GJets.root')
  qcd.AddFile(basedir+'QCD.root')
else:
  data.AddFile(basedir+'MET.root') 
  if args.shower.upper()=='HERWIG':
    prong1.AddFile(basedir+'prong1.root')
    prong2.AddFile(basedir+'prong2_Herwig.root')
    prong3.AddFile(basedir+'prong3_Herwig.root')
    label += 'herwig_'
    print label
  else:
    prong1.AddFile(basedir+'prong1.root')
    prong2.AddFile(basedir+'prong2.root')
    prong3.AddFile(basedir+'prong3.root')
#  prong1.AddFile(basedir+'QCD.root')
#  prong1.AddFile(basedir+'WJets.root')
#  prong2.AddFile(basedir+'Diboson.root')
#  prong2.AddFile(basedir+'TTbar.root')
#  prong2.AddFile(basedir+'SingleTop.root')
#  prong3.AddFile(basedir+'TTbar.root')
#  prong3.AddFile(basedir+'SingleTop.root')

#processes = [prong1,data]
for p in processes:
  plot.AddProcess(p)

plot.AddDistribution(root.Distribution('fj1MSD',50,450,40,'fatjet m_{SD} [GeV]','Events/10 GeV'))
# plot.AddDistribution(root.Distribution('fj1MSDL2L3',50,450,40,'L2L3-corr fatjet m_{SD} [GeV]','Events/10 GeV'))
#plot.AddDistribution(root.Distribution("1",0,2,1,"dummy","dummy"))

### DRAW AND CATALOGUE ###
plot.DrawAll(figsdir+'/'+label)
