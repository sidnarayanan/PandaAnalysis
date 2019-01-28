#!/usr/bin/env python

from PandaCore.Tools.script import *

args = parse('--outdir', '--f', ('--obs', FLOAT))

from PandaCore.Tools.root_interface import Selector
import numpy as np
from re import sub
from glob import glob
from array import array 
load('PandaCoreDrawers')

fpaths = glob(args.f + '/*quant*root')
quantiles = []
for fp in fpaths:
    f = root.TFile.Open(fp)
    t = f.Get('limit')
    t.GetEntry(0)
    fp_ = fp.split('/')[-1]
    fp_ = sub('higgsCombine.*quant','',fp_)
    fp_ = sub('.root','',fp_)
    quantiles.append((float(fp_), t.limit))

quantiles.sort(key = lambda x : x[0])

x = []
y = []
for q,l in quantiles:
    if len(x) == 0 or l > x[-1]:
        x.append(l)
        y.append(q)

y = np.array(y)
x = np.array(x)
dy = np.divide(y[1:]-y[:-1], x[1:]-x[:-1])

x = array('f', x[:-1])
y = array('f', y[:-1])
dy = array('f', dy)

g = root.TGraph(len(x), x, y)
h = root.TH1D('','',20,0,1)
for ib in xrange(1, h.GetNbinsX()+1):
    xval = h.GetBinCenter(ib)
    if xval <= x[0] or xval >= x[:-1]:
        h.SetBinContent(ib, 0)
    else:
        h.SetBinContent(ib, g.Eval(xval))

plot = root.CanvasDrawer()
plot.SetTDRStyle()
plot.InitLegend(.55, .78, .93, .91, 2)
c = plot.GetCanvas()

#h.Scale(1./h.Integral())
#h.Draw('hist')
g.Draw('')

plot.GetLegend().AddEntry(h, '#mu_{0.95} distribution','l')

#l = root.TLine()
#l.SetLineColor(root.kBlack)
#l.SetLineWidth(2)
#l.SetLineStyle(2)
#l.DrawLine(args.obs, 0, args.obs, 0.5)
#
#txt = root.TLatex()
#txt.SetNDC(False)
#txt.SetTextAlign(21)
#txt.DrawLatex(args.obs, 0.7, '#tilde{q}^{obs}_{#mu}')
#txt.SetTextAlign(11)
#txt.DrawLatex(0.2, 5, '#mu=%.2f'%args.mu)
#txt.DrawLatex(0.2, 2.5, 'CL_{s}=%.2f'%cls)
#
#print cls,'=',hmufill.Integral(),'/',h0fill.Integral()

plot.Draw(args.outdir+'/', 'cls_exp')

