#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import linspace
import numpy as np

basedir = getenv('PANDA_FLATDIR') + '/' 

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--proc',metavar='proc',type=str)
parser.add_argument('--outdir',metavar='outdir',type=str)
args = parser.parse_args()

figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Tools.Load import Load
from PandaCore.Tools.Misc import *
from PandaCore.Tools.root_interface import Selector
from math import sqrt
Load('CanvasDrawer')
root.gROOT.SetBatch()
root.gStyle.SetNumberContours(999);
root.gStyle.SetPalette(root.kCool)

SMOOTH = True

plot = root.CanvasDrawer()
plot.SetTDRStyle()
root.gStyle.SetPadRightMargin(0.2)
c = root.TCanvas()
plot.SetCanvas(c)

s = Selector()
s_nlo = Selector()

cut = 'trueGenBosonPt>50 && lheHT>100'
if args.proc == 'WJets':
    cut = tAND(cut, 'genBosonMass>70 && genBosonMass<90')
else:
    cut = tAND(cut, 'genBosonMass>80 && genBosonMass<100')

branches = ['normalizedWeight', '0.001*trueGenBosonPt', '0.001*lheHT']

infile = basedir + args.proc + '.root'
f = root.TFile.Open(infile); t = f.Get('events')
print 'LO'
s.read_tree(t, branches = branches, cut = cut)

infile = basedir + args.proc + '_nlo.root'
f_nlo = root.TFile.Open(infile); t_nlo = f_nlo.Get('events')
print 'NLO'
s_nlo.read_tree(t_nlo, branches = branches, cut = cut)
print 'done reading'

fout = root.TFile.Open('../../data/higgs/hbb_kfactors.root', 'UPDATE')

def draw(x, xbins,  xlabel, y, ybins, ylabel):
    x_ = '0.001*'+x
    y_ = '0.001*'+y
    plot.Reset()
#    plot.AddCMSLabel()
    plot.AddSqrtSLabel()
    xbins = array('f', xbins)
    ybins = array('f', ybins)
    h = root.TH2D('h','h',len(xbins)-1,xbins,len(ybins)-1,ybins)
    hlo = s.draw([x_, y_], weight='normalizedWeight', hbase=h)
    hnlo = s_nlo.draw([x_, y_], weight='normalizedWeight', hbase=h)
    hratio = hnlo.Clone(); hratio.Divide(hlo)
    if SMOOTH:
        hratio.Smooth()
    root.gStyle.SetPaintTextFormat('.2g') 
    hratio.GetXaxis().SetTitle(xlabel)
    hratio.GetYaxis().SetTitle(ylabel)
    hratio.GetZaxis().SetTitle('NLO/LO')
    hratio.Draw('colz text')
    plot.Draw(args.outdir+'/','%s_%s_ratio_%s'%(x, y, args.proc))
    fout.WriteTObject(hratio, 'h_' + args.proc, 'OVERWRITE')
    

draw(
    'trueGenBosonPt', [0.001 * x for x in [50,100,150,200,300,400,500,600,800,1000,1400]], 'p_{T}^{V} [TeV]',
    'lheHT', [0.001 * x for x in [100,200,400,600,1400]], 'H_{T} [TeV]'
    )
  

fout.Close()
