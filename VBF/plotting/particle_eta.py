#!/usr/bin/env python

from PandaCore.Tools.script import *

args = parse('--outdir', '--infile')

from PandaCore.Tools.root_interface import Selector
import numpy as np

root.gROOT.ProcessLine('#include "PandaAnalysis/Utilities/src/PackingHelperStandalone.cc"')

load('PandaTreeObjects')
load('PandaCoreDrawers')
load('PandaAnalysisUtilities')

f = root.TFile.Open(args.infile)
t = f.Get('events')

PFCand = root.panda.PFCand

hbase = root.TH1D('','',20,-5,5)

nmax = root.std.numeric_limits('Short_t').max()
s = Selector()
cut = 'pa::PHS::up(pfCandidates.packedPt)>5'
eta = 'pfCandidates.packedEta * %f / %i'%(6, nmax)
ptype = 'pfCandidates.ptype'
s.read_tree(t, branches=[(eta, -9, 200), (ptype, -100, 200)], cut=cut)
s.flatten()

def get_hist(*ptypes):
    masks = [s.data[ptype] == p for p in ptypes]
    mask = np.zeros_like(s.data[eta]).astype(bool)
    for p in ptypes:
        mask = np.logical_or(mask, s.data[ptype] == p)
    print ptypes
    print s.data[eta][mask].min(), s.data[eta][mask].max()
    h = s.draw(eta, hbase=hbase, mask=mask)
    print h.Integral()
    return h

cands = [
        ('e^{#pm}', [PFCand.ep, PFCand.em]),
        ('#mu^{#pm}', [PFCand.mup, PFCand.mum]),
        ('h^{#pm}', [PFCand.hp, PFCand.hm]),
        ('h^{0}', [PFCand.h0, PFCand.h_HF]),
        ('#gamma', [PFCand.gamma, PFCand.egamma_HF]),
        ]

plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.InitLegend(0.3, 0.65, 0.93, 0.92, 3)
plot.Stack(True)

hists = {}
for i, (label, ptypes) in enumerate(cands):
    h = get_hist(*ptypes)
    h.GetXaxis().SetTitle('Particle #eta')
    h.GetYaxis().SetTitle('Number of particles')
    hists[label] = h

def draw(logy):
    plot.Reset()
    plot.AddSqrtSLabel()
    if logy:
        plot.Logy()
        plot.SetAbsMin(5e3)
    for i, (label, ptypes) in enumerate(cands):
        plot.AddHistogram(hists[label], label, i+5)

draw(False)
plot.Draw(args.outdir+'/', 'particle_eta')

draw(True)
plot.Draw(args.outdir+'/', 'particle_eta_logy')
