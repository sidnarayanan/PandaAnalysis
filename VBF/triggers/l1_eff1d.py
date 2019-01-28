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
parser.add_argument('--finor', action='store_true')
parser.add_argument('--iso', action='store_true')
args = parser.parse_args()

figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Utils.load import Load
from PandaCore.Tools.Misc import *
from PandaCore.Tools.root_interface import Selector
from math import sqrt
Load('GraphAsymmErrDrawer')
root.gROOT.SetBatch()

lumi=36000
import PandaAnalysis.VBF.PandaSelection as sel 

cut = 'fabs(jotEta[{0}])>2.75 && fabs(jotEta[{0}])<3 && filter==1'
cut = tAND(cut, '!(fabs(jotEta[{0}]+2.81)<0.4 && fabs(jotPhi[{0}]-2.07)<0.4)')
if args.finor:
    cut = tAND(cut, 'nJotEC==1')
sigcut = 'finor[1]!=0' if args.finor else 'jotL1EGBX[{0}]==-1'
if args.iso:
    sigcut = tAND(sigcut, 'jotL1EGIso[{0}]>0')

plot = root.GraphAsymmErrDrawer()
plot.SetLineWidth(3)
plot.SetTDRStyle()
plot.InitLegend(.2,.65,.48,.9)

plotlog = root.HistogramDrawer()
plotlog.SetTDRStyle()
plotlog.InitLegend(.6,.65,.88,.9)
plotlog.Logy()

s = Selector()

infile_tmpl = basedir + '%s.root'
labels = ['MET', 'JetHT', 'SingleMuon']
#labels = ['MET']
files = {l:root.TFile.Open(infile_tmpl%l) for l in labels}
trees = {l:files[l].Get('events') for l in labels}

def fn(hbase, branch_tmpl, xtitle, postfix, save=False):
    hnum = {}; hden = {}; hratio = {}
    for label,t in trees.iteritems():
        hnum[label] = hbase.Clone()
        hden[label] = hbase.Clone()
        for iJ in xrange(3):
            branches = [branch_tmpl.format(iJ)]
            s.read_tree(t, branches = branches, cut=tAND(cut, 'jotPt[{0}]>0').format(iJ))
            hden[label].Add(s.draw(branches, hbase=hbase))
            s.read_tree(t, branches = branches, cut=tAND(cut, tAND('jotPt[{0}]>0', sigcut)).format(iJ))
            hnum[label].Add(s.draw(branches, hbase=hbase))

        hratio[label] = root.TEfficiency(hnum[label], hden[label]).CreateGraph()

        hnum[label].GetXaxis().SetTitle(xtitle)
        if args.finor:
            hnum[label].GetYaxis().SetTitle('Number of BX-1 FinOR')
        else:
            hnum[label].GetYaxis().SetTitle('Number of BX-1 matched jets')
        hden[label].GetXaxis().SetTitle(xtitle)
        hden[label].GetYaxis().SetTitle('Number of jets')
        hratio[label].SetMinimum(0)
        hratio[label].SetMaximum(1.5)
        hratio[label].GetXaxis().SetTitle(xtitle)
        hratio[label].SetLineWidth(3)
        if args.finor:
            hratio[label].GetYaxis().SetTitle('#epsilon(BX_{-1} L1A)')
        elif args.iso:
            hratio[label].GetYaxis().SetTitle('IsoEG30 BX=-1 eff')
        else:
            hratio[label].GetYaxis().SetTitle('IsoEG30,EG40 BX=-1 eff')

    suffix = postfix
    if args.finor:
        suffix += '_finor'
    elif args.iso:
        suffix += '_egiso'
    else:
        suffix += '_eg'

    plotlog.Reset()
#    plotlog.AddCMSLabel()
    for i,label in enumerate(trees):
        plotlog.AddHistogram(hnum[label],label,root.kData,i+1,'elp')
    plotlog.Draw(args.outdir+'/','oned_'+suffix+'_num')

    plotlog.Reset()
#    plotlog.AddCMSLabel()
    for i,label in enumerate(trees):
        plotlog.AddHistogram(hden[label],label,root.kData,i+1,'elp')
    plotlog.Draw(args.outdir+'/','oned_'+suffix+'_den')

    plot.Clear()
    plot.ClearLegend()
#    plot.AddCMSLabel()
    for i,label in enumerate(trees):
        plot.AddGraph(hratio[label],label,i+1,1,'lz')
    plot.Draw(args.outdir+'/','oned_'+suffix+'_ratio')

    if save:
        fout = root.TFile.Open('../../data/vbf16/trig/l1.root','update')
        for label,h in hratio.iteritems():
            fout.WriteTObject(h, 'h1_'+label+'_'+suffix, 'overwrite')



hbase = root.TH1D('h0', 'h', 20, 40, 600)
fn(hbase, 'jotPt[{0}]', 'p_{T} [GeV]', 'jotPt', False)

hbase = root.TH1D('h1', 'h', 20, 0, 600)
fn(hbase, 'jotPt[{0}]*jotNEMF[{0}]', 'EM p_{T} [GeV]', 'jotEMPt')

# hbase = root.TH1D('h2', 'h', 20, 30, 256)
# fn(hbase, 'jotL1Pt[{0}]', 'L1 p_{T} [GeV]', 'jotL1Pt')
 
# hbase = root.TH1D('h3', 'h', 20, 2.6, 3.2)
# fn(hbase, 'fabs(jotEta[{0}])', '#eta', 'jotEta')
 
