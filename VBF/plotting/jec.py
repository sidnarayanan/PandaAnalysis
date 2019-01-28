#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange
import numpy as np
from PandaCore.Tools.script import * 
from PandaCore.Tools.root_interface import Selector

args = parse('--infile', '--jet', '--outdir')

figsdir = args.outdir
Load('HistogramDrawer')

plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.SetLumi(36)
#plot.SetAutoRange(False)
plot.InitLegend(.2, .7, .93, .92, 3)

ptarray = array('f', [30, 50, 80, 120, 160, 200, 250, 300, 500, 1000])
etaarray = array('f', [-5, -4, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5])
uncarray = array('f', np.linspace(-30, 30, 200))

f = root.TFile.Open(args.infile)
t = f.Get('events')
ks = t.GetListOfBranches()
pt = args.jet + 'Pt'
rawpt = args.jet + 'RawPt'
eta = 'fabs(' + args.jet + 'Eta)'
corr = '100*(%s-%s)/%s'%(pt, rawpt, rawpt)

s = Selector()
s.read_tree(t, branches=[pt, rawpt, eta, corr], cut='%s>30'%(pt))

def draw(arr, xlabel, xvar, etabins):
  plot.Reset()
  plot.AddLumiLabel()
  plot.GetCanvas().SetLogx('Pt' in xvar)
  plot.cd()

#  plot.Logy()
  hbase = root.TH2D('', '', len(arr)-1, arr, len(uncarray)-1, uncarray)
  h1s = []
  for ie in xrange(len(etabins)-1):
    mask = np.logical_and(s[eta]>etabins[ie], s[eta]<etabins[ie+1])
    h2 = s.draw([xvar, corr], hbase=hbase, mask=mask)
    h1 = h2.QuantilesX(0.5, 'eta_%i'%ie)
    h1.GetXaxis().SetTitle(xlabel)
    h1.GetYaxis().SetTitle('Relative correction [%]')
    plot.AddHistogram(h1, '%.1f < |#eta| < %.1f'%(etabins[ie], etabins[ie+1]), root.kExtra1, ie+1) 

  plot.Draw(args.outdir, 'jec_%s%s'%(args.jet, xvar) )


draw(ptarray, 'Uncorrected jet p_{T} [GeV]', rawpt, [0, 1.5, 2.5, 5])
