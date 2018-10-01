#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
from glob import glob
import argparse
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--var',type=str)
args = parser.parse_args()

figsdir = args.outdir
var=args.var
basedir = getenv('PANDA_FITSCAN')
argv=[]

import ROOT as root
from PandaCore.Utils.load import *

##Color palette

Load('Drawers','GraphAsymmErrDrawer')

plot = root.GraphAsymmErrDrawer()
plot.SetTDRStyle()
plot.SetGrid()
plot.SetLumi(36.6)
plot.AddCMSLabel()
plot.AddLumiLabel(False)
plot.GetCanvas().SetLeftMargin(0.15)

c = root.TCanvas()

cuts = set([])
listoffiles = glob(basedir+'/higgsCombine*.root')
for f in listoffiles:
  try:
    fname = f.split('/')[-1]
    fname = fname.replace('higgsCombine','').replace('.Asymptotic.mH120.root','')
    cut = float(fname.split('_')[1])
    cuts.add(cut)
  except IndexError:
    pass

def xformSet(s):
  s_ = sorted(list(s))
  return s_
  l = [(3*s_[0]-s_[1])/2]
  l += [(s_[x]+s_[x+1])/2 for x in xrange(len(s_)-1)]
  l += [(3*s_[-1]-s_[-2])/2]
  return l

cutarray = array('f',xformSet(cuts))
limitarray = array('f',[0 for x in cutarray])
limitdownarray = array('f',[0 for x in cutarray])
limituparray = array('f',[0 for x in cutarray])
zeros = array('f',[0 for x in cutarray])

for iC in xrange(len(cuts)):
    try:
      cut = cutarray[iC]
      fname = basedir+'/higgsCombine_%.2f_.Asymptotic.mH120.root'%(cut)
      f = root.TFile(fname)
      t = f.Get('limit')
      t.GetEntry(2)
      limitarray[iC] = t.limit
      t.GetEntry(1)
      limitdownarray[iC] = -t.limit+limitarray[iC]
      t.GetEntry(3)
      limituparray[iC] = t.limit-limitarray[iC]
    except:
      pass

c.Clear(); c.cd()
g = root.TGraphAsymmErrors(len(cuts),cutarray,limitarray,zeros,zeros,zeros,zeros)
gerr = root.TGraphAsymmErrors(len(cuts),cutarray,limitarray,zeros,zeros,limitdownarray,limituparray)
gerr.GetXaxis().SetTitle('BDT cut')
gerr.GetYaxis().SetTitle('#sigma_{95% CL}/#sigma_{theory}')
gerr.GetYaxis().SetTitleOffset(1.35)
gerr.SetMinimum(0); g.SetMaximum(0.15)
plot.SetLineWidth(3)
g.SetFillColor(root.kGreen+1)
gerr.SetLineWidth(0)
gerr.SetFillColorAlpha(root.kGreen+1,0.5)
plot.AddGraph(gerr,'',0,1,'3')
plot.AddGraph(g,'',0,1,'l2')
plot.Draw(args.outdir+'/','optimized_scan_'+var)

