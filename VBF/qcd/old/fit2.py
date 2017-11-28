#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange, sqrt

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default='/home/snarayan/public_html/figs/vbf/v13/qcd_estimate/')
args = parser.parse_args()
infile = args.outdir+'qcd_met_mjj.root'

figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import Load 
Load('HistogramDrawer')

f = root.TFile.Open(infile)
h = f.Get('h_mjj_met')

double_exp = root.TF2('double_exp', '([0]*TMath::Exp(-[3]*y))*TMath::Exp(-([1]*TMath::Exp(-[4]*y))*x) + ([2]*TMath::Exp(-[5]*y))*TMath::Exp(-[3]*x)', 50, 1000)
# double_exp.SetParNames('A_{1}', 'k_{1}', 'A_{2}', 'k_{2}', 'q_{1}', 'p_{1}')
NPARAMS=6
double_exp.SetLineColor(root.kRed)

init = {
    0 : 2500,
    1 : 0.003,
    2 : 50,
    3 : 0.07,
    # 4 : 0,
}

for k,v in init.iteritems():
    double_exp.SetParameter(k,v)

h_params = [h.ProjectionX('h_%i'%iP) for iP in xrange(NPARAMS)]
for hh in h_params:
    hh.Clear()

root.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(int(1e9))
plot = root.HistogramDrawer()
plot.SetTDRStyle()
# root.gStyle.SetPadRightMargin(0.15)
c = root.TCanvas()
c.cd() 
c.SetLogy()
plot.Logy(True)
plot.InitLegend()
plot.SetCanvas(c)
root.gStyle.SetOptStat(0)

xaxis = h.GetXaxis()
yaxis = h.GetYaxis()

for iMjj in xrange(1, h.GetNbinsX()):
    xwidth = xaxis.GetBinWidth(iMjj)
    for iMet in xrange(1, h.GetNbinsY()):
        ywidth = yaxis.GetBinWidth(iMet)
        area = xwidth * ywidth
        h.SetBinContent(iMjj, iMet, h.GetBinContent(iMjj, iMet) / area)
        h.SetBinError(iMjj, iMet, h.GetBinError(iMjj, iMet) / area)

for _ in xrange(3):
    h.Fit(double_exp, 'ME')


h_projections = []
for iMjj in xrange(1, h.GetNbinsX()+1):
    h_proj = h.ProjectionY('metProjection_%i'%iMjj, iMjj, iMjj)

    plot.Reset(True)
    plot.AddCMSLabel()
    plot.AddSqrtSLabel()
    xaxis = h.GetXaxis()
    lo = xaxis.GetBinLowEdge(iMjj)
    hi = lo + xaxis.GetBinWidth(iMjj)
    plot.AddPlotLabel('%.0f < m_{jj} < %.0f [GeV]'%(lo, hi), .18, .8, False, 42, .04)
    plot.AddHistogram(h_proj, 'Monte Carlo', root.kData)

    fit_proj = root.TF12('f12_%i'%iMjj, double_exp, (lo+hi)/2, 'x')
    plot.AddAdditional(fit_proj, 'l', 'Fit')

    
    plot.Draw(args.outdir, 'fit2d_met_%i'%iMjj)

for iMet in xrange(1, h.GetNbinsY()+1):
    h_proj = h.ProjectionX('mjjProjection_%i'%iMet, iMet, iMet)

    plot.Reset(True)
    plot.AddCMSLabel()
    plot.AddSqrtSLabel()
    xaxis = h.GetXaxis()
    lo = xaxis.GetBinLowEdge(iMet)
    hi = lo + xaxis.GetBinWidth(iMet)
    plot.AddPlotLabel('%.0f < E_{T}^{miss} < %.0f [GeV]'%(lo, hi), .18, .8, False, 42, .04)
    plot.AddHistogram(h_proj, 'Monte Carlo', root.kData)

    fit_proj = root.TF12('f12_mjj_%i'%iMet, double_exp, (lo+hi)/2, 'y')
    plot.AddAdditional(fit_proj, 'l', 'Fit')

    
    plot.Draw(args.outdir, 'fit2d_mjj_%i'%iMet)