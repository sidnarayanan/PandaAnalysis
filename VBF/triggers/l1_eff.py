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
parser.add_argument('--spike',action='store_true')
parser.add_argument('--finor',action='store_true')
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

cut = 'fabs(jotEta[{0}])>2.4 || (jotNEMF[{0}]<0.9 && jotNHF[{0}]<0.8) && filter==1 && jotPt[{0}]>0'
if args.vbf:
    cut = tAND(cut, tAND('nJot>1', sel.mjj))
elif args.spike:
    cut = tAND(cut, '!(fabs(jotL1Eta[{0}]+2.81)<0.2 && fabs(jotL1Phi[{0}]-2.07)<0.2)')
if args.finor:
    cut = tAND(cut, 'nJotEC==1 && fabs(jotEta[{0}])>2.25 && fabs(jotEta[{0}])<3')
sigcut = 'finor[1]!=0' if args.finor else 'jotL1EGBX[{0}]==-1'

plot = root.CanvasDrawer()
plot.SetTDRStyle()
plot.AddPlotLabel(args.f.replace('_',' ').replace('AllEras','Run2016*'), .18, .77, False, 42, .04)
if args.vbf:
    plot.AddPlotLabel('#Delta#phi_{jj}<1.5, #Delta#eta_{jj}>1', .18, .72, False, 42, .04)
elif args.spike:
    plot.AddPlotLabel('|#phi-2.01|>0.2 or |#eta+2.81|>0.2', .18, .72, False, 42, .04)
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
        s.read_tree(t, branches = branches, cut=cut.format(iJ))
        hden.Add(s.draw(branches, hbase=hbase))
        s.read_tree(t, branches = branches, cut=tAND(cut, sigcut).format(iJ))
        hnum.Add(s.draw(branches, hbase=hbase))

    hratio = hnum.Clone()
    eff = root.TEfficiency(hnum, hden)
    for ix in xrange(1, hratio.GetNbinsX()+1):
        for iy in xrange(1, hratio.GetNbinsY()+1):
            hratio.SetBinContent(ix, iy, eff.GetEfficiency(hratio.GetBin(ix, iy)))

    if args.vbf:
        suffix = 'vbf'
    elif args.spike:
        suffix = 'spike'
    else:
        suffix = 'inclusive'
    if args.finor:
        suffix += '_finor'
    suffix += postfix

    plot.Reset()
    plot.AddCMSLabel()
    hnum.GetXaxis().SetTitle(xtitle)
    hnum.GetYaxis().SetTitle(ytitle)
    hnum.GetZaxis().SetTitle('Number of BX-1 jets')
    hnum.Draw('colz')
    plot.Draw(args.outdir+'/',args.f+'_'+suffix+'_num')

    plot.Reset()
    plot.AddCMSLabel()
    hden.GetXaxis().SetTitle(xtitle)
    hden.GetYaxis().SetTitle(ytitle)
    hden.GetZaxis().SetTitle('Number of jets')
    hden.Draw('colz')
    plot.Draw(args.outdir+'/',args.f+'_'+suffix+'_den')

    plot.Reset()
    plot.AddCMSLabel()
    hratio.SetMinimum(0)
    hratio.SetMaximum(1)
    hratio.GetXaxis().SetTitle(xtitle)
    hratio.GetYaxis().SetTitle(ytitle)
    hratio.GetZaxis().SetTitle('L1IsoEG BX=-1 eff')
    hratio.Draw('colz')
    plot.Draw(args.outdir+'/',args.f+'_'+suffix)
    if save:
        fout = root.TFile.Open('../../data/vbf16/trig/l1.root','update')
        fout.WriteTObject(hratio, 'h_'+args.f+'_'+suffix, 'overwrite')
        fout.Close()
        f.cd()

hbase = root.TH2D('h3', 'h', 20, -3.1416, 3.1416, 20, -5, 5)
fn(hbase, ['jotPhi[{0}]', 'jotEta[{0}]'], '#phi', '#eta', '_etaphi', False)

hbase = root.TH2D('h-2', 'h', 20, 0, 1000, 20, 30, 600)
fn(hbase, ['met', 'fabs(jotPt[{0}])'], 'p_{T}^{miss} [GeV]', 'p_{T} [GeV]', '_metpt', False)

hbase = root.TH2D('h-1', 'h', 20, 0, 1000, 20, 0, 5)
fn(hbase, ['met', 'fabs(jotEta[{0}])'], 'p_{T}^{miss} [GeV]', '|#eta|', '_meteta', False)

hbase = root.TH2D('h0', 'h', 20, 40, 600, 20, 0, 5)
fn(hbase, ['jotPt[{0}]', 'fabs(jotEta[{0}])'], 'p_{T} [GeV]', '|#eta|', '', True)

hbase = root.TH2D('h1', 'h', 20, 40, 600, 20, 0, 5)
fn(hbase, ['jotNEMF[{0}]*jotPt[{0}]', 'fabs(jotEta[{0}])'], 'p_{T}^{EM} [GeV]', '|#eta|', '_emf', True)

hbase = root.TH2D('h8', 'h', 20, 0, 1, 20, 0, 5)
fn(hbase, ['jotNEMF[{0}]', 'fabs(jotEta[{0}])'], 'EM Fraction', '|#eta|', '_emfrac', False)
fn(hbase, ['jotNHF[{0}]', 'fabs(jotEta[{0}])'], 'Hadron Fraction', '|#eta|', '_hfrac', False)

hbase = root.TH2D('h2', 'h', 20, 0, 300, 20, 0, 5)
fn(hbase, ['jotL1Pt[{0}]', 'fabs(jotL1Eta[{0}])'], 'p_{T}^{L1} [GeV]', '|#eta_{L1}|', '_l1', False)

