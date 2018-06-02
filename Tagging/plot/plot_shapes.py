#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
import argparse

basedir = getenv('PANDA_FLATDIR')
figsdir = basedir+'/figs'

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--indir',metavar='indir',type=str,default=basedir)
parser.add_argument('--outdir',metavar='outdir',type=str,default=figsdir)
parser.add_argument('--cut',metavar='cut',type=str,default=None)
args = parser.parse_args()

figsdir = args.outdir
basedir = args.indir
argv=[]

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
tcut = root.TCut
Load('Drawers','PlotUtility')

### SET GLOBAL VARIABLES ###
lumi = 12918.
logy=False
cut = 'nFatJet==1 && fj1Pt>250 && nLooseLep==1 && UWmag>250' # so we select leptonic wjets or semileptonic tt

if not args.cut:
  label = 'noCut_'
  plotlabel = '40 GeV < m_{SD}'
  cut = tAND(cut,'fj1MSD>40')
elif args.cut=='mass':
  cut = tAND(cut,'fj1MSD>110 && fj1MSD<210')
  label = 'massCut_'
  plotlabel = '110 < m_{SD} < 210 GeV'
elif args.cut=='massW':
  cut = tAND(cut,'fj1MSD>50 && fj1MSD<100')
  label = 'massWCut_'
  plotlabel = '50 < m_{SD} < 100 GeV'

### LOAD PLOTTING UTILITY ###
plot = root.PlotUtility()
plot.Stack(False)
plot.InitLegend()
plot.SetCut(tcut(cut))
plot.SetLumi(lumi/1000)
plot.SetNormFactor(True)
plot.SetTDRStyle()
plot.AddCMSLabel()
# plot.AddLumiLabel()
if plotlabel:
  plot.AddPlotLabel(plotlabel,.18,.77,False,42,.04)

weight = '%f*normalizedWeight*sf_pu*sf_lep*sf_ewkV*sf_qcdV*sf_tt'%lumi
plot.SetMCWeight(weight)

### DEFINE PROCESSES ###
wjetsq     = root.Process('W+q',root.kExtra7); wjetsq.additionalCut = root.TCut('abs(fj1HighestPtGen)!=21')
wjetsg     = root.Process('W+g',root.kExtra2); wjetsg.additionalCut = root.TCut('abs(fj1HighestPtGen)==21')
ttbar     = root.Process('t#bar{t} [matched]',root.kTTbar); ttbar.additionalCut = root.TCut('(fj1IsMatched==1&&fj1GenSize<1.44)')
ttbarunmatched     = root.Process('t#bar{t} [unmatched]',root.kExtra1); ttbarunmatched.additionalCut = root.TCut('(fj1IsMatched==0||fj1GenSize>1.44)')
processes = [wjetsg,wjetsq,ttbarunmatched,ttbar]

### ASSIGN FILES TO PROCESSES ###
wjetsq.AddFile(basedir+'WJets.root')
wjetsg.AddFile(basedir+'WJets.root')
ttbar.AddFile(basedir+'TTbar.root')
ttbarunmatched.AddFile(basedir+'TTbar.root')

for p in processes:
  plot.AddProcess(p)

plot.AddDistribution(root.Distribution('top_ecfv7_bdt',-1.2,1.,20,'Top ECF BDT v7','Events'))

plot.AddDistribution(root.Distribution('fj1MSD',40,450,20,'fatjet m_{SD} [GeV]','Events/12.5 GeV'))
#plot.AddDistribution(root.Distribution('top_ecf_bdt',-0.5,.5,20,'Top ECF BDT','Events'))

plot.AddDistribution(root.Distribution('fj1Tau32SD',0,1,20,'Groomed #tau_{32}','Events',999,-999,'tau32SD'))

plot.AddDistribution(root.Distribution('fj1Tau32',0,1,20,'#tau_{32}','Events',999,-999,'tau32'))

plot.AddDistribution(root.Distribution('jet1Pt',15,500,20,'leading jet p_{T} [GeV]','Events'))

plot.AddDistribution(root.Distribution('nJet',-0.5,8.5,9,'N_{jet}','Events'))

plot.AddDistribution(root.Distribution('UWmag',250,500,20,'W recoil [GeV]','Events'))

plot.AddDistribution(root.Distribution('puppimet',0,750,20,'MET [GeV]','Events/37.5 GeV'))

plot.AddDistribution(root.Distribution('fj1Pt',250,1000,20,'fatjet p_{T} [GeV]','Events/37.5 GeV'))


### DRAW AND CATALOGUE ###
plot.DrawAll(figsdir+'/'+label)
