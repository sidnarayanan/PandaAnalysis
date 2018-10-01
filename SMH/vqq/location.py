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
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
parser.add_argument('--prefix',type=str,default='inclusive_')
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
c.SetLogz()
plot.SetCanvas(c)

s = Selector()
f = root.TFile.Open(basedir + '/' + 'VqqHbb' + '.root')
t = f.Get('events')
weight = sel.weights['signal']%40000

branches = [weight]
for b in ['Eta', 'Phi']:
    for i,j in [('fjHiggsIdx', 'H'), ('fjVIdx', 'V')]:
        branches.append('fjQ[%s]'%i)
        branches.append('fj%s[%s]'%(b, i))
        branches.append('gen%s%s'%(j, b))
        branches.append('fj%s[%s] - gen%s%s'%(b, i, j, b))

s.read_tree(t, branches=branches, cut=tAND(sel.cuts['signal'], args.cut))

hbase = root.TH2D('h','h', 20, -2.5, 2.5, 20, -3.14159, 3.14159)
def draw(fmt, label_fmt, fname):
    plot.Reset()
    h = s.draw([fmt.format('Eta'), fmt.format('Phi')], weight=weight, hbase=hbase)
    h.GetXaxis().SetTitle(label_fmt.format('#eta'))
    h.GetYaxis().SetTitle(label_fmt.format('#phi'))
    h.Draw('colz')
    plot.Draw(args.outdir+'/',args.prefix+fname)
    
draw('fj{0}[fjHiggsIdx]', 'Reco H {0}', 'hReco')
draw('fj{0}[fjVIdx]', 'Reco V {0}', 'vReco')
draw('genH{0}', 'Gen H {0}', 'hGen')
draw('genV{0}', 'Gen V {0}', 'vGen')

hbase = root.TH2D('h2','h', 20, 0, 2.5, 20, 0, 2*3.14159)
draw('fj{0}[fjHiggsIdx] - genH{0}', '|Reco H {0} - Gen H {0}|', 'hDiff')
draw('fj{0}[fjVIdx] - genV{0}', '|Reco V {0} - Gen V {0}|', 'vDiff')

hbase = root.TH2D('h3','h', 20, -1, 1, 20, -1, 1)
plot.Reset()
h = s.draw(['fjQ[fjHiggsIdx]', 'fjQ[fjVIdx]'], weight=weight, hbase=hbase)
h.GetXaxis().SetTitle('Reco H Q')
h.GetYaxis().SetTitle('Reco V Q')
h.Draw('colz')
plot.Draw(args.outdir+'/',args.prefix+'charge')
