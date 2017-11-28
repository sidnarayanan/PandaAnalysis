#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange, sqrt

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default='/home/snarayan/public_html/figs/vbf/v13/qcd_estimate2/')
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

# double_exp = root.TF1('double_exp', "([0]-[2]/[3]*(TMath::Exp(-[3]*[6])-TMath::Exp(-[3]*[7]))-[4]/[5]*(TMath::Exp(-[5]*[6])-TMath::Exp(-[5]*[7])))*([1]/(TMath::Exp(-[1]*[6])-TMath::Exp(-[1]*[7])))*TMath::Exp(-[1]*x) + [2]*TMath::Exp(-[3]*x) + [4]*TMath::Exp(-[5]*x)", 50, 1000)
# double_exp.SetParNames("Integral", "k_{1}", "A_{2}", "k_{2}", "A_{3}", "k_{3}", "E_{1}", "E_{2}")
double_exp = root.TF1('double_exp', "[0]*TMath::Exp(-[1]*x) + [2]*TMath::Exp(-[3]*x)", 50, 1000)
double_exp.SetParNames('A_{1}', 'k_{1}', 'A_{2}', 'k_{2}')
double_exp.SetLineColor(root.kRed)
NPARAM = double_exp.GetNpar()

h_params = [h.ProjectionX('h_%i'%iP) for iP in xrange(NPARAM)]
for hh in h_params:
    hh.Clear()

root.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(int(1e9))
plot = root.HistogramDrawer()
plot.SetTDRStyle()
c = root.TCanvas()
c.cd() 
c.SetLogy()
plot.Logy(True)
plot.InitLegend()
plot.SetCanvas(c)
root.gStyle.SetOptStat(0)

parameters = {x:[] for x in xrange(NPARAM)}
paramerrors = {x:[] for x in xrange(NPARAM)}

h_projections = []
for iMjj in xrange(1, h.GetNbinsX()+1):
    h_proj = h.ProjectionY('metProjection_%i'%iMjj, iMjj, iMjj)
    integral = h_proj.Integral()
    for iMet in xrange(1, h_proj.GetNbinsX()+1):
        h_proj.SetBinContent(iMet, h_proj.GetBinContent(iMet)/h_proj.GetBinWidth(iMet))
        h_proj.SetBinError(iMet, h_proj.GetBinError(iMet)/h_proj.GetBinWidth(iMet))
        h_proj.SetBinError(iMet, sqrt(pow(h_proj.GetBinError(iMet), 2) + pow(0.1*h_proj.GetBinContent(iMet), 2)))
    double_exp.SetParameter(0, integral)
    double_exp.SetParameter(1, 0.1)
    double_exp.SetParameter(2, .1*integral)
    double_exp.SetParameter(3, 0.01)
    double_exp.SetParLimits(0, 0.1*integral, 10*pow(integral, 3))
    double_exp.SetParLimits(1, 0, 100)
    double_exp.SetParLimits(2, 0, integral)
    double_exp.SetParLimits(3, 0, 10)
#     A2_estimation = root.TMath.Exp(0.02*300.)*h_proj.GetBinContent(h_proj.FindBin(300.));
#     A3_estimation = root.TMath.Exp(0.01*700.)*h_proj.GetBinContent(h_proj.FindBin(700.));
#     double_exp.FixParameter(0, integral);
#     double_exp.SetParameter(1, 0.05);
#     double_exp.SetParLimits(1, 0.00001, 0.5);
#     double_exp.SetParameter(2, A2_estimation );
#     double_exp.SetParLimits(2, 0.01*A2_estimation, 100.*A2_estimation);
#     double_exp.SetParameter(3, 0.02);
#     double_exp.SetParLimits(3, 0.00001, 0.5);
#     double_exp.SetParameter(4, A3_estimation );
#     double_exp.SetParLimits(4, 0, A2_estimation);
#     double_exp.SetParameter(5, 0.005);
#     double_exp.SetParLimits(5, 0.00001, 0.1);
#     double_exp.FixParameter(6, 50);
#     double_exp.FixParameter(7, 1000);
# #    double_exp.FixParameter(0, integral)
# #    double_exp.SetParameter(1, 0.05)
# #    double_exp.SetParameter(2, 0.001*integral)
# #    double_exp.SetParameter(3, 1)
# #    double_exp.SetParameter(4, 1)
# #    for p in [3,4]:
# #        double_exp.SetParLimits(p, 0, 20)
#     # double_exp.FixParameter(4, 0.)
#     # double_exp.SetParameter(5, 1)
#     # double_exp.SetParLimits(0, integral, 10*pow(integral, 3))
#     # double_exp.SetParLimits(1, 0.001, 0.2)
#     # double_exp.SetParLimits(2, 0.01, 0.95*integral)
#     # double_exp.SetParLimits(3, 0, 10)
#     # double_exp.SetParLimits(4, 0., 0.)
#     # double_exp.SetParLimits(5, -2, 2)

    plot.Reset(True)
    plot.AddCMSLabel()
    plot.AddSqrtSLabel()
    xaxis = h.GetXaxis()
    lo = xaxis.GetBinLowEdge(iMjj)
    hi = lo + xaxis.GetBinWidth(iMjj)
    plot.AddPlotLabel('%.0f < m_{jj} < %.0f [GeV]'%(lo, hi), .18, .8, False, 42, .04)
    plot.AddHistogram(h_proj, 'Monte Carlo', root.kData)

    print 'fitting',iMjj
    for _ in xrange(3):
        h_proj.Fit(double_exp, 'ME', '', 50, 1000)
    plot.AddAdditional(double_exp, 'l', 'Fit')
    for x in xrange(NPARAM):
        parameters[x].append(double_exp.GetParameter(x))
        paramerrors[x].append(double_exp.GetParError(x))

    plot.Draw(args.outdir, 'fit_%i'%iMjj)

plot.Logy(False)
for x in xrange(NPARAM):
    plot.Reset(True)
    plot.AddCMSLabel()
    plot.AddSqrtSLabel()
    h_p = h_params[x]
    for iMjj in xrange(0, h.GetNbinsX()):
        h_p.SetBinContent(iMjj+1, parameters[x][iMjj])
        h_p.SetBinError(iMjj+1, paramerrors[x][iMjj])
    plot.AddHistogram(h_p, 'Fit parameter '+double_exp.GetParName(x), root.kSignal1, root.kSignal1, 'elp')
    plot.Draw(args.outdir,'fit_param%i'%x)
