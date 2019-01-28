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
uncarray = array('f', np.linspace(0, 0.1, 200))

f = root.TFile.Open(args.infile)
t = f.Get('events')
ks = t.GetListOfBranches()
pt = args.jet + 'Pt'
pt_ = pt + '_'
uncs = []
for k in ks:
  if k.GetName().startswith(pt_) and k.GetName().endswith('Up'):
    uncs.append(k.GetName().replace(pt_, '')[:-2])

s = Selector()
branches = {u: 'fabs(%s-%s%sUp)/%s'%(pt,pt_,u,pt) for u in uncs}
all_branches = [pt, args.jet+'Eta'] + branches.values()
s.read_tree(t, branches=all_branches, cut='%s>30'%(pt))

c2 = root.TCanvas()

majors = ['Absolute', 'Relative', 'PileUp', 'SinglePion', 'Flavor', 'Fragmentation', 'Time', 'Total']

def draw(arr, xlabel, xvar):
  plot.Reset()
  plot.AddLumiLabel()
  plot.GetCanvas().SetLogx('Pt' in xvar)
#  plot.Logy()
  hbase = root.TH2D('', '', len(arr)-1, arr, len(uncarray)-1, uncarray)
  h2s = {u : s.draw([xvar, branches[u]], hbase=hbase) for u in uncs}
  h1s = {u : h.QuantilesX(0.68, u+'_qx') for u,h in h2s.iteritems()}
  norms = {u : h.Integral() for u,h in h1s.iteritems() if 'Total' not in u}
  norms = sorted(norms, key=lambda x : norms[x])
  norms = norms[-7:]
  maxy = 1.6*max([h.GetMaximum() for _,h in h1s.iteritems()])

  plot.cd()

  hms = {m:h1s.values()[0].Clone(m+'_clone') for m in majors}
  for _,h in hms.iteritems():
    for ib in xrange(1, h.GetNbinsX()+1):
      h.SetBinContent(ib, 0)

  htotal = hms['Total']
  for u in uncs:
    if 'Total' in u:
      continue
    h = h1s[u]
    for m in majors:
      if m in u:
        break
    else:
      print 'could not classify',u
      m = 'Other'
    hm = hms[m]
    for ib in xrange(1, h.GetNbinsX()+1):
      hm.SetBinContent(ib, hm.GetBinContent(ib) + (h.GetBinContent(ib)**2))
      htotal.SetBinContent(ib, htotal.GetBinContent(ib) + (h.GetBinContent(ib)**2))

  for i,u in enumerate(majors):
    hm.GetXaxis().SetTitle(xlabel)
    hm.GetYaxis().SetTitle('Relative uncertainty [%]')
#    hm.SetMaximum(maxy)
    hm = hms[u]
    for ib in xrange(1, h.GetNbinsX()+1):
      hm.SetBinContent(ib,  100 * hm.GetBinContent(ib)**0.5)
    if 'Total' in u:
      utotal = u
      continue
    plot.AddHistogram(hm, u, 5, i+2)

  
  hms[utotal].SetFillColorAlpha(root.kGray, 0.5)
  hms[utotal].SetLineColorAlpha(root.kBlack, 0.5)
  hms[utotal].SetLineStyle(3)
  plot.AddHistogram(hms[utotal], utotal, 5, root.kBlack)

  plot.Draw(args.outdir, 'var_%s%s'%(args.jet, xvar) )


draw(ptarray, 'Jet p_{T} [GeV]', args.jet+'Pt')
draw(etaarray, '#eta', args.jet+'Eta')
