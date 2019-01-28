#!/usr/bin/env python

from PandaCore.Tools.script import *

args = parse('--outdir', '--f', ('--mu', FLOAT), ('--obs', FLOAT))

from PandaCore.Tools.root_interface import Selector
import numpy as np
load('PandaCoreDrawers')

f = root.TFile.Open(args.f)
t = f.Get('q') 

plot = root.CanvasDrawer()
plot.SetTDRStyle()
plot.InitLegend(.55, .78, .93, .91, 2)
c = plot.GetCanvas()
c.SetLogy()

hbase = root.TH1D('','',50,0,10)

s = Selector()
def get_hist(typ, lc):
    s.read_tree(t, branches=['q'], cut='type==%i'%typ)
    h = s.draw('q', hbase=hbase)
    h.GetXaxis().SetTitle('#tilde{q}_{#mu}')
    h.GetYaxis().SetTitle('p(#tilde{q}_{#mu}|#bar{#mu})')
    h.SetLineColor(lc)
    h.SetLineWidth(2)
    h.SetMinimum(1)
    return h

hmu = get_hist(1, root.kRed) 
h0 = get_hist(-1, root.kBlue)

hmu.Scale(1./hmu.Integral())
h0.Scale(1./h0.Integral())
hmu.SetMaximum(10)

hmufill = hmu.Clone('mufill')
h0fill = h0.Clone('0fill')
for ib in xrange(1,hbase.GetNbinsX()+1):
    if hbase.GetBinCenter(ib) < args.obs:
        hmufill.SetBinContent(ib, 0)
        h0fill.SetBinContent(ib, 0)
h0fill.SetFillColorAlpha(root.kBlue,0.5)
hmufill.SetFillColorAlpha(root.kRed,0.5)
hmufill.SetLineWidth(0)
h0fill.SetLineWidth(0)

cls = hmufill.Integral() / h0fill.Integral()

hmu.Draw('hist')
h0.Draw('hist same')
h0fill.Draw('hist same')
hmufill.Draw('hist same')

plot.GetLegend().AddEntry(hmu, '#bar{#mu}=#mu','l')
plot.GetLegend().AddEntry(h0, '#bar{#mu}=0', 'l')
plot.GetLegend().AddEntry(hmufill, 'CL_{s+b}','f')
plot.GetLegend().AddEntry(h0fill, 'CL_{b}', 'f')

l = root.TLine()
l.SetLineColor(root.kBlack)
l.SetLineWidth(2)
l.SetLineStyle(2)
l.DrawLine(args.obs, 0, args.obs, 0.5)

txt = root.TLatex()
txt.SetNDC(False)
txt.SetTextAlign(21)
txt.DrawLatex(args.obs, 0.7, '#tilde{q}^{obs}_{#mu}')
txt.SetTextAlign(11)
txt.DrawLatex(0.2, 5, '#mu=%.2f'%args.mu)
txt.DrawLatex(0.2, 2.5, 'CL_{s}=%.2f'%cls)

print cls,'=',hmufill.Integral(),'/',h0fill.Integral()

plot.Draw(args.outdir+'/', 'cls_%.2f'%args.mu)

