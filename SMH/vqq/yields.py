#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import linspace
import numpy as np

basedir = getenv('PANDA_FLATDIR') + '/' 

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--cut',metavar='cut',type=float)
args = parser.parse_args()

figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Utils.load import Load
from PandaCore.Tools.Misc import *
from PandaCore.Tools.root_interface import Selector
import PandaAnalysis.SMH.VqqSel as sel 
from math import sqrt
Load('CanvasDrawer')
root.gROOT.SetBatch()
root.gStyle.SetNumberContours(999);
root.gStyle.SetPalette(root.kCool)

SMOOTH = True

plot = root.CanvasDrawer()
plot.SetTDRStyle()
plot.SetLumi(40000)
root.gStyle.SetPadRightMargin(0.2)
c = root.TCanvas()
plot.SetCanvas(c)

selectors = {}
files = {}
trees = {}
for p in ['QCD','VqqHbb']:
    selectors[p] = {}
    files[p] = root.TFile.Open(basedir + '/' + p + '.root')
    trees[p] = files[p].Get('events')
    for cat in ['pass', 'fail']:
        selectors[p][cat] = Selector()
weight = sel.weights['signal']%40000

branches = [weight, 'fjMSD[fjHiggsIdx]', 'fjMSD[fjVIdx]']

for p,cats in selectors.iteritems():
    for cat,s in cats.iteritems():
        logger.info('sensitivity.py', '%s %s'%(p, cat))
        print sel.cuts[cat]%(args.cut, args.cut)
        s.read_tree(trees[p], branches=branches, cut=sel.cuts[cat]%(args.cut, args.cut))

hbase = root.TH2D('h','h', 10, 30, 300, 10, 30, 300)
def draw(cat):
    plot.Reset()
#    plot.AddSqrtSLabel()
    h = {p : selectors[p][cat].draw(['fjMSD[fjHiggsIdx]', 'fjMSD[fjVIdx]'], 
                                    weight=weight, hbase=hbase)
         for p in ['VqqHbb', 'QCD']}
    hsens = hbase.Clone()
    for iX in xrange(1, hsens.GetNbinsX()+1):
        for iY in xrange(1, hsens.GetNbinsY()+1):
            s = h['VqqHbb'].GetBinContent(iX, iY)
            b = h['QCD'].GetBinContent(iX, iY)
            hsens.SetBinContent(iX, iY, s / (np.sqrt(b) + 0.1 * b))
    root.gStyle.SetPaintTextFormat('.2g') 
    hsens.GetXaxis().SetTitle('m_{H} [GeV]')
    hsens.GetYaxis().SetTitle('m_{V} [GeV]')
    hsens.GetZaxis().SetTitle('S/(#sqrt{B}+0.1B)')
    hsens.Draw('colz text')
    plot.Draw(args.outdir+'/','%s_%.1f'%(cat, args.cut))
    

draw('pass')
draw('fail')
