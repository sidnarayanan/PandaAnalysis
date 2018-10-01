#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_ZEYNEPDIR')+'/merged/' 
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cut',metavar='cut',type=str,default='inc')
args = parser.parse_args()
sname=argv[0]
argv=[]
lumi = 12900.

import ROOT as root
from ROOT import gROOT
from PandaCore.Utils.load import *
from PandaCore.Tools.Misc import *
from PandaAnalysis.VBF.LooseSelection import *
from array import array
from math import sqrt
Load('Drawers','HistogramDrawer')

### DEFINE REGIONS ###
vbfsel = noid
if args.cut=='vbf':
  vbfsel = tAND(noid,'jjDEta>4.5')
if args.cut=='vbf_loose':
  vbfsel = tAND(noid,'jjDEta>3.5')
label = args.cut

recoilBins = [200,250,300,350,400,500,600,1000]
nRecoilBins = len(recoilBins)-1
recoilBins = array('f',recoilBins)

logy=False
plot = root.HistogramDrawer()
plot.SetLumi(12.9)
plot.SetTDRStyle()
plot.AddCMSLabel()
plot.AddLumiLabel()
plot.SetNormFactor(True)
plot.Logy(logy)
plot.SetAbsMin(0.0001)
plot.SetTDRStyle()
plot.InitLegend()

plotr = root.HistogramDrawer()
plotr.SetLumi(12.9)
plotr.SetRatioStyle()
plotr.AddCMSLabel()
plotr.AddLumiLabel()
plotr.InitLegend(.15,.6,.5,.8)

counter=0

# fvbf = root.TFile(baseDir+'monojet_VBF_HToInvisible_M125_13TeV_powheg_pythia8.root'); tvbf = fvbf.Get('events')
fqcdz = root.TFile(baseDir+'ZJets.root'); tqcdz = fqcdz.Get('events')
fewkz = root.TFile(baseDir+'EWKZJets.root'); tewkz = fewkz.Get('events')
fqcdw = root.TFile(baseDir+'WJets.root'); tqcdw = fqcdw.Get('events')
fewkw = root.TFile(baseDir+'EWKWJets.root'); tewkw = fewkw.Get('events')
fqcdzvv = root.TFile(baseDir+'ZtoNuNu.root'); tqcdzvv = fqcdzvv.Get('events')
fewkzvv = root.TFile(baseDir+'EWKZtoNuNu.root'); tewkzvv = fewkzvv.Get('events')

ctmp = root.TCanvas()

def getMETHist(tree,cut,binwidth=False):
  global counter
  ctmp.cd()
  h = root.TH1D('h%i'%counter,'h%i'%counter,nRecoilBins,recoilBins)
  h.GetXaxis().SetTitle('U')
  h.GetYaxis().SetTitle('arbitrary units')
  if binwidth:
    h.Scale(1,'width')
  tree.Draw('met>>h%i'%counter,cut)
  counter += 1
  h.SetFillStyle(0)
  return h

def plotRatio(V,label):
  if V=='Z':
    tqcdsig = tqcdzvv
    tewksig = tewkzvv
    tqcdbg = tqcdz
    tewkbg = tewkz
  else:
    tqcdsig = tqcdw
    tewksig = tewkw
    tqcdbg = tqcdw
    tewkbg = tewkw
  hsigqcd = getMETHist(tqcdsig,tTIMES(tAND(vbfsel,cuts['signal']),weights['signal']))
  hsigewk = getMETHist(tewksig,tTIMES(tAND(vbfsel,cuts['signal']),weights['signal']))
  hsigewk.Divide(hsigqcd)
  if V=='Z':
    hmuqcd = getMETHist(tqcdbg,tTIMES(tAND(vbfsel,cuts['zmm']),weights['zmm']))
    hmuewk = getMETHist(tewkbg,tTIMES(tAND(vbfsel,cuts['zmm']),weights['zmm']))
    heleqcd = getMETHist(tqcdbg,tTIMES(tAND(vbfsel,cuts['zee']),weights['zee']))
    heleewk = getMETHist(tewkbg,tTIMES(tAND(vbfsel,cuts['zee']),weights['zee']))
  else:
    hmuqcd = getMETHist(tqcdbg,tTIMES(tAND(vbfsel,cuts['wmn']),weights['wmn']))
    hmuewk = getMETHist(tewkbg,tTIMES(tAND(vbfsel,cuts['wmn']),weights['wmn']))
    heleqcd = getMETHist(tqcdbg,tTIMES(tAND(vbfsel,cuts['wen']),weights['wen']))
    heleewk = getMETHist(tewkbg,tTIMES(tAND(vbfsel,cuts['wen']),weights['wen']))
  hmuewk.Add(heleewk)
  hmuqcd.Add(heleqcd)
  hmuewk.Divide(hmuqcd)
  for h in [hsigewk,heleewk,hmuewk]:
    h.GetYaxis().SetTitle('EWK %s/QCD %s'%(V,V))
  hsigerr = hsigewk.Clone()
  hsigerr.SetFillStyle(3004)
  hsigerr.SetFillColorAlpha(root.kBlack,0.5)
  hsigerr.SetLineWidth(0)
  hmuerr = hmuewk.Clone()
  hmuerr.SetFillStyle(3005)
  hmuerr.SetFillColorAlpha(root.kMagenta+1,0.5)
  hmuerr.SetLineWidth(0)
  plotr.AddHistogram(hsigewk,'Signal region',root.kExtra2,root.kBlack)
  if V=='W':
    plotr.AddHistogram(hmuewk,'Single lep CRs',root.kExtra1)
  #  plotr.AddHistogram(heleewk,'Single e CR',root.kExtra4)
  else:
    plotr.AddHistogram(hmuewk,'Dilep CRs',root.kExtra1)
  #  plotr.AddHistogram(heleewk,'Dielectron CR',root.kExtra4)
  plotr.AddAdditional(hsigerr,'e2')
  plotr.AddAdditional(hmuerr,'e2')
  plotr.Draw(args.outdir,label)
  plotr.Reset()

def plotWZ(label):
  totalcut = tTIMES(tAND(vbfsel,cuts['signal']),weights['signal'])
  hqcdz = getMETHist(tqcdzvv,totalcut)
  hqcdw = getMETHist(tqcdw,totalcut)
  hewkz = getMETHist(tewkzvv,totalcut)
  hewkw = getMETHist(tewkw,totalcut)
  hewkw.Divide(hqcdw)
  hewkz.Divide(hqcdz)
  for h in [hewkw,hewkz]:
    h.GetYaxis().SetTitle('EWK V/QCD V')
  hwerr = hewkw.Clone()
  hwerr.SetFillStyle(3004)
  hwerr.SetFillColorAlpha(root.kGreen+1,0.5)
  hwerr.SetLineWidth(0)
  hzerr = hewkz.Clone()
  hzerr.SetFillStyle(3005)
  hzerr.SetFillColorAlpha(root.kCyan+1,0.5)
  hzerr.SetLineWidth(0)
  plotr.AddHistogram(hewkz,'EWK Z/QCD Z',root.kExtra2)
  plotr.AddHistogram(hewkw,'EWK W/QCD W',root.kExtra3)
  plotr.AddAdditional(hzerr,'e2')
  plotr.AddAdditional(hwerr,'e2')
  plotr.Draw(args.outdir,label)
  plotr.Reset()


plotWZ('%s_WZ_ratio'%label)
plotRatio('Z','%s_Z_ratio'%label)
plotRatio('W','%s_W_ratio'%label)
