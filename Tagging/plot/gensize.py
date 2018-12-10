#!/usr/bin/env python

from PandaCore.Tools.script import * 
from re import sub 
from collections import namedtuple
import numpy as np
from sys import exit 
from PandaCore.Tools.root_interface import Selector

Load('CanvasDrawer')
c = root.CanvasDrawer(800, 600)
c.SetTDRStyle()
c.GetCanvas().SetRightMargin(0.125)
c.GetCanvas().SetTopMargin(0.05)
c.GetCanvas().SetLeftMargin(0.13)
c.GetCanvas().SetBottomMargin(0.15)
c.cd()


args = parse('--outdir')

f = root.TFile.Open(flatdir + '/Top.root')
t = f.Get('events')
s = Selector()

s.read_tree(t, branches=['ptweight','gen_pt','TMath::Sqrt(gen_size)'],
            cut='fabs(gen_eta)<2.5')

max_x = 500
max_y = 3.14159
hbase = root.TH2D('','',50,50,max_x,50,0,max_y)
h = s.draw(fields=['gen_pt','TMath::Sqrt(gen_size)'],
           hbase=hbase, weight='ptweight')
h.GetXaxis().SetTitle('Top p_{T} [GeV]')
h.GetYaxis().SetTitle("max#DeltaR_{qq}")
h.GetXaxis().SetTitleOffset(1); h.GetXaxis().SetTitleSize(0.06)
h.GetYaxis().SetTitleOffset(1); h.GetYaxis().SetTitleSize(0.06)

def get_q(q, c, y):
    hq = h.QuantilesX(q)
    hq.SetLineColor(c)
    hq.SetLineWidth(3)
    intersection = None
    return hq, hq.GetBinContent(hq.FindBin(250))
    for ib in xrange(hq.GetNbinsX()+1,0,-1):
        if hq.GetBinContent(ib) > y:
            intersection = hq.GetBinCenter(ib)
            break
    return hq, intersection

hq5, int5 = get_q(0.5, root.kRed, 1.2)
hq7, int7 = get_q(0.75, root.kBlack, 1.5)

#h2.SetMinimum(0)
#h2.SetMinimum(-15)
h.Draw('colz')
hq5.Draw('hist same')
#hq7.Draw('hist same')

def draw_line(intersection, y):
    l = root.TLine()
    l.SetLineColor(root.kGray)
    l.SetLineWidth(2)
    l.SetLineStyle(2)
    # l.DrawLine(0,y,intersection,y)
    # l.DrawLine(intersection,0,intersection,y)
    l.DrawLine(50,intersection,250,intersection)
    l.DrawLine(250,0,250,intersection)

draw_line(int5, 250)
#draw_line(int5, 1.2)
#draw_line(int7, 1.5)

c.Draw(args.outdir, 'ptdr')

