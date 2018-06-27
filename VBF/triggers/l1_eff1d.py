#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange
import numpy as np

basedir = getenv('PANDA_FLATDIR') + '/' 

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
args = parser.parse_args()

figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Utils.load import Load
from PandaCore.Tools.Misc import *
from PandaCore.Tools.root_interface import Selector
from math import sqrt
Load('HistogramDrawer')
root.gROOT.SetBatch()

lumi=36000
import PandaAnalysis.VBF.PandaSelection as sel 

cut = 'fabs(jotEta[{0}])>2.75 && fabs(jotEta[{0}])<3 && jotPt[{0}]>100'
cut = tAND(cut, '!(fabs(jotEta[{0}]+2.75)<0.25 && fabs(jotPhi[{0}]-2.1)<0.25)')

plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.InitLegend()

plotlog = root.HistogramDrawer()
plotlog.SetTDRStyle()
plotlog.InitLegend()
plotlog.Logy()

s = Selector()

infile_tmpl = basedir + '%s_AllEras.root'
labels = ['MET', 'JetHT', 'SingleMuon']
files = {l:root.TFile.Open(infile_tmpl%l) for l in labels}
trees = {l:files[l].Get('events') for l in labels}

def fn(hbase, branch_tmpl, xtitle, postfix):
    hnum = {}; hden = {}; hratio = {}
    for label,t in trees.iteritems():
        hnum[label] = hbase.Clone()
        hden[label] = hbase.Clone()
        for iJ in xrange(3):
            branches = [branch_tmpl.format(iJ)]
            s.read_tree(t, branches = branches, cut=tAND(cut, 'jotPt[{0}]>0').format(iJ))
            hden[label].Add(s.draw(branches, hbase=hbase))
            s.read_tree(t, branches = branches, cut=tAND(cut, 'jotL1EGBX[{0}]==-1').format(iJ))
            hnum[label].Add(s.draw(branches, hbase=hbase))

        hratio[label] = hnum[label].Clone()
        for ib in xrange(1, hnum[label].GetNbinsX()+1):
            den = hden[label].GetBinContent(ib)
            if den > 0:
                hratio[label].SetBinContent(ib, hnum[label].GetBinContent(ib) / den)
                hratio[label].SetBinError(ib, hnum[label].GetBinError(ib) / den)

        hnum[label].GetXaxis().SetTitle(xtitle)
        hnum[label].GetYaxis().SetTitle('Number of BX-1 matched jets')
        hden[label].GetXaxis().SetTitle(xtitle)
        hden[label].GetYaxis().SetTitle('Number of jets')
        hratio[label].SetMinimum(0)
        hratio[label].SetMaximum(1)
        hratio[label].GetXaxis().SetTitle(xtitle)
        hratio[label].GetYaxis().SetTitle('L1IsoEG BX=-1 eff')

    suffix = postfix

    plotlog.Reset()
    plotlog.AddCMSLabel()
    for i,label in enumerate(trees):
        plotlog.AddHistogram(hnum[label],'%s [%i evts]'%(label, hnum[label].Integral()),root.kData,i+1,'elp')
    plotlog.Draw(args.outdir+'/','oned_'+suffix+'_num')

    plotlog.Reset()
    plotlog.AddCMSLabel()
    for i,label in enumerate(trees):
        plotlog.AddHistogram(hden[label],'%s [%i evts]'%(label, hden[label].Integral()),root.kData,i+1,'elp')
    plotlog.Draw(args.outdir+'/','oned_'+suffix+'_den')

    plot.Reset()
    plot.AddCMSLabel()
    for i,label in enumerate(trees):
        plot.AddHistogram(hratio[label],label,root.kData,i+1,'elp')
    plot.Draw(args.outdir+'/','oned_'+suffix+'_ratio')


hbase = root.TH1D('h0', 'h', 30, 40, 800)
fn(hbase, 'jotPt[{0}]', 'p_{T} [GeV]', 'jotPt')

hbase = root.TH1D('h1', 'h', 30, 0, 800)
fn(hbase, 'jotPt[{0}]*jotNEMF[{0}]', 'EM p_{T} [GeV]', 'jotEMPt')

hbase = root.TH1D('h2', 'h', 20, 30, 256)
fn(hbase, 'jotL1Pt[{0}]', 'L1 p_{T} [GeV]', 'jotL1Pt')

hbase = root.TH1D('h3', 'h', 20, 2.6, 3.2)
fn(hbase, 'fabs(jotEta[{0}])', '#eta', 'jotEta')

