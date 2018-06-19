#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange
import numpy as np

basedir = getenv('PANDA_FLATDIR') + '/' 

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--f',metavar='f',type=str)
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--vbf',action='store_true')
args = parser.parse_args()

figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Utils.load import Load
from PandaCore.Tools.Misc import *
from PandaCore.Tools.root_interface import Selector
from math import sqrt
Load('CanvasDrawer')
root.gROOT.SetBatch()
root.gStyle.SetNumberContours(999);
root.gStyle.SetPalette(root.kCool)

lumi=36000
import PandaAnalysis.VBF.PandaSelection as sel 

cut = tAND('nJot>1 && met>250', sel.mjj) if args.vbf else '1==1'

plot = root.CanvasDrawer()
plot.SetTDRStyle()
plot.AddPlotLabel(args.f.replace('_',' ').replace('AllEras','Run2016*'), .18, .77, False, 42, .04)
if args.vbf:
    plot.AddPlotLabel('#Delta#phi_{jj}<1.5, #Delta#eta_{jj}>1', .18, .72, False, 42, .04)
root.gStyle.SetPadRightMargin(0.2)
c = root.TCanvas()
plot.SetCanvas(c)

s = Selector()

infile = basedir + args.f + '.root'
f = root.TFile.Open(infile); t = f.Get('events')

def fn(hbase, branch_tmpl, xtitle, ytitle, postfix, save=False):
    hnum = hbase.Clone()
    hden = hbase.Clone()
    for iJ in xrange(3):
        branches = [x.format(iJ) for x in branch_tmpl]
        print branches
        s.read_tree(t, branches = branches, cut=tAND(cut, 'jotPt[%i]>0'%(iJ)))
        hden.Add(s.draw(branches, hbase=hbase))
        s.read_tree(t, branches = branches, cut=tAND(cut, 'jotL1EGBX[%i]==-1'%(iJ)))
        hnum.Add(s.draw(branches, hbase=hbase))

    hratio = hnum.Clone()
    hratio.Divide(hden)

    plot.Reset()
    plot.AddCMSLabel()
    hratio.SetMinimum(0)
    hratio.SetMaximum(1)
    hratio.GetXaxis().SetTitle(xtitle)
    hratio.GetYaxis().SetTitle(ytitle)
    hratio.GetZaxis().SetTitle('L1IsoEG BX=-1 eff')
    hratio.Draw('colz')
    suffix = 'vbf' if args.vbf else 'inclusive'
    suffix += postfix
    plot.Draw(args.outdir+'/',args.f+'_'+suffix)
    if save:
        fout = root.TFile.Open('../../data/vbf16/trig/l1.root','update')
        fout.WriteTObject(hratio, 'h_'+args.f+'_'+suffix)
        fout.Close()
        f.cd()

hbase = root.TH2D('h0', 'h', 20, 40, 600, 20, 0, 5)
fn(hbase, ['jotPt[{0}]', 'fabs(jotEta[{0}])'], 'p_{T} [GeV]', '|#eta|', '', True)

hbase = root.TH2D('h1', 'h', 20, 40, 600, 20, 0, 5)
fn(hbase, ['jotNEMF[{0}]*jotPt[{0}]', 'fabs(jotEta[{0}])'], 'p_{T}^{EM} [GeV]', '|#eta|', '_emf', False)

hbase = root.TH2D('h2', 'h', 20, 0, 300, 20, 0, 5)
fn(hbase, ['jotL1Pt[{0}]', 'fabs(jotEta[{0}])'], 'p_{T}^{TP} [GeV]', '|#eta_{TP}|', '_l1', False)

hbase = root.TH2D('h3', 'h', 20, -3.1416, 3.1416, 20, -5, 5)
fn(hbase, ['jotPhi[{0}]', 'jotEta[{0}]'], '#phi', '#eta', '_etaphi', False)

