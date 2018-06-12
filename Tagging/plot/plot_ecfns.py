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
from PandaCore.Utils.load import *
from PandaCore.Tools.Misc import *
tcut = root.TCut
Load('Drawers','PlotUtility')

### SET GLOBAL VARIABLES ###
lumi = 12918.
logy=False
cut = 'nFatjet==1 && fj1Pt>250 && nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && UWmag>250 && fj1MSD>40'
weight = '%f*normalizedWeight*sf_pu*sf_lep*sf_ewkV*sf_qcdV*sf_tt'%lumi
plotlabel=None
if args.cut=='tag':
  cut = tAND(cut,'fj1MaxCSV>0.46 && isojetNBtags==1')
  weight = tTIMES(weight,'sf_sjbtag1*sf_btag1')
  label = 'tag_'
elif args.cut=='mistag':
  cut = tAND(cut,'fj1MaxCSV<0.46 && isojetNBtags==0')
  weight = tTIMES(weight,'sf_sjbtag0*sf_btag0')
  label = 'mistag_'

### LOAD PLOTTING UTILITY ###
plot = root.PlotUtility()
plot.Stack(True)
plot.InitLegend()
plot.SetCut(tcut(cut))
plot.Ratio(1) 
plot.FixRatio(.5)
plot.SetLumi(lumi/1000)
plot.DrawMCErrors(True)
plot.SetTDRStyle()
plot.SetNormFactor(False)
plot.AddCMSLabel()
plot.AddLumiLabel()
if plotlabel:
  plot.AddPlotLabel(plotlabel,.18,.77,False,42,.04)

plot.SetMCWeight(weight)

### DEFINE PROCESSES ###
wjetsq     = root.Process('W+q',root.kWjets); wjetsq.additionalCut = root.TCut('abs(fj1HighestPtGen)!=21')
wjetsg     = root.Process('W+g',root.kExtra2); wjetsg.additionalCut = root.TCut('abs(fj1HighestPtGen)==21')
diboson   = root.Process('Diboson',root.kDiboson)
ttbar     = root.Process('t#bar{t} [matched]',root.kTTbar); ttbar.additionalCut = root.TCut('(fj1IsMatched==1&&fj1GenSize<1.44)')
ttbarunmatched     = root.Process('t#bar{t} [unmatched]',root.kExtra1); ttbarunmatched.additionalCut = root.TCut('(fj1IsMatched==0||fj1GenSize>1.44)')
singletop = root.Process('Single t',root.kST)
qcd       = root.Process("QCD",root.kQCD)
data      = root.Process("Data",root.kData); data.additionalCut = root.TCut('(trigger&1)!=0')
processes = [diboson,singletop,wjetsg,wjetsq,ttbarunmatched,ttbar]

### ASSIGN FILES TO PROCESSES ###
wjetsq.AddFile(basedir+'WJets.root')
wjetsg.AddFile(basedir+'WJets.root')
diboson.AddFile(basedir+'Diboson.root')
ttbar.AddFile(basedir+'TTbar.root')
ttbarunmatched.AddFile(basedir+'TTbar.root')
singletop.AddFile(basedir+'SingleTop.root')
qcd.AddFile(basedir+'QCD.root')

data.AddFile(basedir+'MET.root') 
processes.append(data)

for p in processes:
  plot.AddProcess(p)

plot.AddDistribution(root.Distribution('fj1MSD',40,450,20,'fatjet m_{SD} [GeV]','Events'))

ecfns = [
    ('fj1ECFN_1_2_05',0,.5),
    ('fj1ECFN_1_2_10',0,.6),
    ('fj1ECFN_1_2_20',0,.8),
    ('fj1ECFN_1_2_40',0,1),
    ('fj1ECFN_1_3_05',0,.12),
    ('fj1ECFN_1_3_10',0,.1),
    ('fj1ECFN_1_3_20',0,.1),
    ('fj1ECFN_1_3_40',0,.15),
    ('fj1ECFN_2_3_05',0,.15),
    ('fj1ECFN_2_3_10',0,.15),
    ('fj1ECFN_2_3_20',0,.25),
    ('fj1ECFN_2_3_40',0,.2),
    ('fj1ECFN_3_3_05',0,.2),
    ('fj1ECFN_3_3_10',0,.25),
    ('fj1ECFN_3_3_20',0,.3),
    ('fj1ECFN_3_3_40',0,.3),
    ('fj1ECFN_1_4_05',0,.05),
    ('fj1ECFN_1_4_10',0,.05),
    ('fj1ECFN_1_4_20',0,.05),
    ('fj1ECFN_1_4_40',0,.05),
    ('fj1ECFN_2_4_05',0,.02),
    ('fj1ECFN_2_4_10',0,.05),
    ('fj1ECFN_2_4_20',0,.05),
    ('fj1ECFN_2_4_40',0,.05),
  ]

def convert(f):
  f_ = f.replace('fj1ECFN_','').split('_')
  return 'e(%s,%s,%.1f)'%(f_[0],f_[1],float(f_[2])/10.)

for ecfn in ecfns:
  var = ecfn[0]
  lo = ecfn[1]
  hi = ecfn[2]
  plot.AddDistribution(root.Distribution(var,lo,hi,20,convert(var),'Events'))

### DRAW AND CATALOGUE ###
plot.DrawAll(figsdir+'/'+label)
