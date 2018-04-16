#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
from glob import glob
import argparse

basedir = getenv('PANDA_VBFSCAN')
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--indir',metavar='indir',type=str,default=basedir)
#parser.add_argument('--outname',type=str)
parser.add_argument('--var1',type=str)
parser.add_argument('--var2',type=str) 
parser.add_argument('--var3',type=str,default=None) # this one is not scanned, just fixed 
args = parser.parse_args()
figsdir = args.outdir
var1=args.var1; var2=args.var2; var3=args.var3
basedir=args.indir
argv=[]

titles = {
    'deta':'min #Delta#eta(j_{1},j_{2})',
    'fabsdphi':'max #Delta#phi(j_{1},j_{2})',
    'mjj':'min M(j_{1},j_{2}) [GeV]',
    }

import ROOT as root
from PandaCore.Tools.Load import *

##Color palette
'''
ncontours = 999;
stops = [     0.10,     0.2000,      0.40,   0.600,     0.800,    1.0000]
red   = [255./255.,  169./255.,  91./255., 36./255.,   8./255.,   6./255.]
green = [255./255.,  194./255., 149./255., 99./255.,  46./255.,   0./255.]
blue  = [255./255.,  243./255., 240./255.,240./255., 240./255., 233./255.]
stopsArray = array('d',stops)
redArray   = array('d',red)
greenArray = array('d',green)
blueArray  = array('d',blue)
root.TColor.CreateGradientColorTable(len(stops), stopsArray, redArray, greenArray, blueArray, ncontours);
'''
root.gStyle.SetNumberContours(999);
root.gStyle.SetPalette(root.kBird)
root.gStyle.SetPaintTextFormat(".3g")


Load('CanvasDrawer')

plot = root.CanvasDrawer()
plot.SetTDRStyle()
root.gStyle.SetPadRightMargin(0.2)

c = root.TCanvas()

cut1s = set([]); cut2s = set([])
if var3:
  listoffiles = glob(basedir+'/higgsCombine%s*%s*%s*.root'%(var1,var2,var3))
else:
  listoffiles = glob(basedir+'/higgsCombine%s*%s*.root'%(var1,var2))
for f in listoffiles:
  try:
    fname = f.split('/')[-1]
    fname = fname.replace('higgsCombine','').replace('.Asymptotic.mH120.root','')
    cut1 = float(fname.split('_')[1])
    cut2 = float(fname.split('_')[4])
    cut1s.add(cut1)
    cut2s.add(cut2)
  except IndexError as e:
    pass

def xformSet(s):
  s_ = sorted(list(s))
  l = [(3*s_[0]-s_[1])/2]
  l += [(s_[x]+s_[x+1])/2 for x in xrange(len(s_)-1)]
  l += [(3*s_[-1]-s_[-2])/2]
  return l

cut1array = array('f',xformSet(cut1s))
cut2array = array('f',xformSet(cut2s))
h2 = root.TH2F('h2','h2',len(cut1array)-1,cut1array,len(cut2array)-1,cut2array)
h2x = h2.GetXaxis(); h2y = h2.GetYaxis();
h2x.SetTitle(titles[var1])
h2y.SetTitle(titles[var2])
h2.GetZaxis().SetTitle('Toy limit')
#h2.SetMaximum(0.5); h2.SetMinimum(0.15)
h2.SetMarkerColor(0)

for e in cut1s:
  for p in cut2s:
    try:
      fname = basedir+'/higgsCombine'
      if var1=='mjj':
        fname += '%s_%i__%s_%.2f'%(var1,int(e),var2,p)
      elif var2=='mjj':
        fname += '%s_%.2f__%s_%i'%(var1,e,var2,int(p))
      else:
        fname += '%s_%.2f__%s_%.2f'%(var1,e,var2,p)
      if var3:
        fname += '__%s'%(var3)
      fname += '.Asymptotic.mH120.root'
      f = root.TFile(fname)
      t = f.Get('limit')
      t.GetEntry(2)
      h2.SetBinContent(h2x.FindBin(e),h2y.FindBin(p),t.limit)
    except:
      pass

c.Clear(); c.cd()
h2.Draw('colz text')
plot.AddCMSLabel(.16,.94)
plot.SetLumi(36.6); plot.AddLumiLabel(True)
plot.SetCanvas(c)
if var3:
  plot.Draw(args.outdir+'/','optimized_scan_'+'_'.join([var1,var2,var3]))
else:
  plot.Draw(args.outdir+'/','optimized_scan_'+'_'.join([var1,var2]))

