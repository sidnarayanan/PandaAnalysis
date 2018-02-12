#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange
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

lumi=36000
import PandaAnalysis.VBF.PandaSelection as sel 

plot = root.CanvasDrawer()
plot.SetTDRStyle()
root.gStyle.SetPadRightMargin(0.2)
c = root.TCanvas()
plot.SetCanvas(c)

s = Selector()

cut = sel.cuts['singlemuon']

infile = basedir + args.proc + '.root'
f = root.TFile.Open(infile); t = f.Get('events')
branches = ['jot2Eta','looseLep1Eta','jot1Eta']
s.read_tree(t, branches = branches, cut = cut)

def draw(x, xbins,  xlabel, y, ybins, ylabel):
    plot.Reset()
    plot.AddCMSLabel()
    plot.AddSqrtSLabel()
    xbins = array('f', xbins)
    ybins = array('f', ybins)
    h = root.TH2D('h','h',len(xbins)-1,xbins,len(ybins)-1,ybins)
    h = s.draw([x, y], weight=None, hbase=h)
    root.gStyle.SetPaintTextFormat('.2g') 
    h.GetXaxis().SetTitle(xlabel)
    h.GetYaxis().SetTitle(ylabel)
    h.GetZaxis().SetTitle('Events')
    h.Draw('colz')
    plot.Draw(args.outdir+'/','compare_%s_%s_%s'%(x, y, args.proc))
    
  
draw(
    'looseLep1Eta', np.arange(-5,5.5,0.5), 'Lepton #eta',
    'jot2Eta', np.arange(-5,5.5,0.5), 'Jet 2 #eta',
    )
 
draw(
    'looseLep1Eta', np.arange(-5,5.5,0.5), 'Lepton #eta',
    'jot1Eta', np.arange(-5,5.5,0.5), 'Jet 1 #eta',
    )
 
draw(
    'jot1Eta', np.arange(-5,5.5,0.5), 'Jet 1 #eta',
    'jot2Eta', np.arange(-5,5.5,0.5), 'Jet 2 #eta',
    )
 
