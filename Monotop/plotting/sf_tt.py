#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
args = parser.parse_args()
args.outdir += '/'
lumi = 36560.
blind=True
linear=False
sname = argv[0]

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Utils.load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
Load('Drawers','HistogramDrawer')


plot = root.HistogramDrawer()
plot.Logy(True)
plot.SetTDRStyle()
plot.InitLegend()
plot.AddCMSLabel()
plot.SetLumi(lumi/1000)
plot.AddLumiLabel(True)
plot.SetNormFactor(True)

plotr = root.HistogramDrawer()
plotr.SetRatioStyle()
plotr.InitLegend()
plotr.AddCMSLabel()
plotr.SetLumi(lumi/1000)
plotr.AddLumiLabel(True)

fin = root.TFile('/home/snarayan/home000/store/panda/v_8026_0_0_gen/TTbar.root')
tree = fin.Get('events')

def build_hist(name,color,width,style=1):
  hnew = root.TH1D(name,name,50,0,1000)
  hnew.GetXaxis().SetTitle('(p_{T}^{t}+p_{T}^{#bar{t}})/2 [GeV]')
  hnew.GetYaxis().SetTitle('(d#sigma/dp_{T})/#sigma')
  hnew.SetLineColor(color)
  hnew.SetLineWidth(width)
  hnew.SetLineStyle(style)
  return hnew

plotstring = 'min(999,max(0,(genTopPt+genAntiTopPt)/2)) >> %s'

integrals = {}
hists = {}

name = 'tt8TeV'
print name
hists[name] = build_hist(name,root.kBlue,3,2)
tree.Draw(plotstring%name,'normalizedWeight*sf_tt8TeV')
integrals[name] = hists[name].Integral()
plot.AddHistogram(hists[name],name,root.kQCD+4,root.kBlue)

name = 'tt8TeV_ext'
print name
hists[name] = build_hist(name,root.kRed,3,2)
tree.Draw(plotstring%name,'normalizedWeight*sf_tt8TeV_ext')
integrals[name] = hists[name].Integral()
plot.AddHistogram(hists[name],name,root.kQCD+5,root.kRed)

'''
name = 'tt8TeV_bound'
print name
hists[name] = build_hist(name,root.kViolet,3,2)
tree.Draw(plotstring%name,'normalizedWeight*sf_tt8TeV_bound')
integrals[name] = hists[name].Integral()
plot.AddHistogram(hists[name],name,root.kQCD+6,root.kViolet)
'''

name = 'tt'
print name
hists[name] = build_hist(name,root.kBlue,3)
tree.Draw(plotstring%name,'normalizedWeight*sf_tt')
integrals[name] = hists[name].Integral()
plot.AddHistogram(hists[name],name,root.kQCD+1,root.kBlue)

name = 'tt_ext'
print name
hists[name] = build_hist(name,root.kRed,3)
tree.Draw(plotstring%name,'normalizedWeight*sf_tt_ext')
integrals[name] = hists[name].Integral()
plot.AddHistogram(hists[name],name,root.kQCD+2,root.kRed)

'''
name = 'tt_bound'
print name
hists[name] = build_hist(name,root.kViolet,3)
tree.Draw(plotstring%name,'normalizedWeight*sf_tt_bound')
integrals[name] = hists[name].Integral()
plot.AddHistogram(hists[name],name,root.kQCD+3,root.kViolet)
'''

name = 'nominal'
print name
hists[name] = build_hist(name,1,3)
tree.Draw(plotstring%name,'normalizedWeight')
integrals[name] = hists[name].Integral()
plot.AddHistogram(hists[name],name,root.kQCD,1)

def print_row(name):
  format_string = '%20s | %10f | %.4f'
  nom_int = integrals['nominal']
  this_int = integrals[name]
  weight = nom_int/this_int
  print format_string%(name,this_int,weight)

print '%20s | %10s | %s'%('Weight','Integral','Normalization')
print '-'*50
print_row('nominal')
print_row('tt')
print_row('tt_ext')
#print_row('tt_bound')
print_row('tt8TeV')
print_row('tt8TeV_ext')
#print_row('tt8TeV_bound')
print '-'*50

plot.Draw(args.outdir,'/top_pt_shapes')

for k,v in hists.iteritems():
  if k=='nominal':
    continue
  v.Divide(hists['nominal'])
  v.SetMaximum(1.5); v.SetMinimum(0)
  v.GetYaxis().SetTitle('#sigma(weighted)/#sigma(nominal)')
hists['nominal'].Divide(hists['nominal'])
hists['nominal'].SetMaximum(1.5)
hists['nominal'].SetMinimum(0)
hists['nominal'].GetYaxis().SetTitle('#sigma(weighted)/#sigma(nominal)')
plotr.AddHistogram(hists['nominal'],'',root.kQCD,1)
plotr.AddHistogram(hists['tt'],'tt',root.kQCD+1,root.kRed)
plotr.AddHistogram(hists['tt_ext'],'tt_ext',root.kQCD+2,root.kBlue)
#plotr.AddHistogram(hists['tt_bound'],'tt',root.kQCD+3,root.kViolet)
plotr.AddHistogram(hists['tt8TeV'],'tt8TeV',root.kQCD+4,root.kRed)
plotr.AddHistogram(hists['tt8TeV_ext'],'tt8TeV_ext',root.kQCD+5,root.kBlue)
#plotr.AddHistogram(hists['tt_bound'],'tt',root.kQCD+6,root.kViolet)

plotr.Draw(args.outdir,'/top_pt_ratios')
