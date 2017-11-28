#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange, sqrt

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,
                    default='/home/snarayan/public_html/figs/vbf/v13/qcd_estimate2/')
args = parser.parse_args()
infile = args.outdir+'qcd_met_mjj.root'

figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import Load 
Load('HistogramDrawer')
root.gInterpreter.GenerateDictionary("std::pair<const std::string, TH1*>", "utility;string;TH1.h")
root.gInterpreter.GenerateDictionary("std::map<std::string, TH1*>", "map;string;TH1.h")
from ROOT import std
def set_item(self, key, value):
    p = std.pair('const std::string, TH1*')(key, value)
    self.insert(p)


f = root.TFile.Open(infile)
h2 = f.Get('h_mjj_met')
NBINSX = h2.GetNbinsX()
#NBINSX = 5
NBINSY = h2.GetNbinsY()

# global fit stuff
met = root.RooRealVar('met','E_{T}^{miss} [GeV]',50,1000)
sample = root.RooCategory('sample', '')
for x in xrange(NBINSX):
    sample.defineType('mjjBin%i'%x)
simult = root.RooSimultaneous('simult','simult',sample)

# build data hist model 
h_data = []
dh_data = []
for x in xrange(NBINSX):
    h1 = h2.ProjectionY('metProj%i'%x, x+1, x+1) 
    for y in xrange(NBINSY):
        h1.SetBinContent(y, h1.GetBinContent(y)/h1.GetBinWidth(y))
        h1.SetBinError(y, h1.GetBinError(y)/h1.GetBinWidth(y))
        if h1.GetBinContent(y) == 0:
            h1.SetBinError(y, 1)
    h_data.append(h1)
    dh_data.append(root.RooDataHist('data%i'%x, 'data%i'%x, root.RooArgList(met), h1))

hmap = root.std.map('std::string, TH1*')()
for x in xrange(NBINSX):
    print h_data[x]
    print type(h_data[x])
    set_item(hmap, 'mjjBin%i'%x, h_data[x])
data  = root.RooDataHist('data', 
                        'data', 
                        root.RooArgList(met),
                        sample,
                        hmap)


## build the fit models
k1 = root.RooRealVar('k1', 'k1', -.09, -2, 0)
k2 = root.RooRealVar('k2', 'k2', -.01, -.02, -0.0)
k3 = root.RooRealVar('k3', 'k3', -0.003, -2, 0)
# k1 = root.RooRealVar('k1', 'k1', -1)
# k2 = root.RooRealVar('k2', 'k2', -1)
# k3 = root.RooRealVar('k3', 'k3', -0.1)
bullshit = [] # prevent python garbage collection from kicking in
def build_3exp(Ns, name):
#     k1 = root.RooRealVar(name+'_k1', name+'_k1', -.095, -2, 0)
#     k2 = root.RooRealVar(name+'_k2', name+'_k2', -.01, -.03, -0.0)
#     k3 = root.RooRealVar(name+'_k3', name+'_k3', -0.003, -2, 0)
    exp1 = root.RooExponential(name+'_exp1', name+'_exp1', met, k1)
    exp2 = root.RooExponential(name+'_exp2', name+'_exp2', met, k2)
    exp3 = root.RooExponential(name+'_exp3', name+'_exp3', met, k3)
    bullshit.extend([exp1, exp2, exp3, k1, k2, k3])
    # model = root.RooAddPdf(name+'_model', name+'_model',
    #                        root.RooArgList(exp1, exp2),
    #                        root.RooArgList(*Ns[:2]))
    model = root.RooAddPdf(name+'_model', name+'_model',
                           root.RooArgList(exp1, exp2, exp3),
                           root.RooArgList(*Ns))
    # if name=='mjjBinModel0':
    #     model.Print('v')
    # return exp1 
    return model 
models = []
params = [k1,k2,k3]
for x in xrange(NBINSX):
    h = h_data[x]
    integral = sum([h.GetBinContent(y)*h.GetBinWidth(y) for y in xrange(1,NBINSY+1)])
    print integral
    N1 = root.RooRealVar('N1%i'%x, 'N1%i'%x, 3*integral, 0.1, 500*integral)
    N2 = root.RooRealVar('N2%i'%x, 'N2%i'%x, 5e-10*integral, 0, 1e-5*integral)
    N3 = root.RooRealVar('N3%i'%x, 'N3%i'%x, 1e-11*integral, 0, 1e-10*integral)
    # N1 = root.RooRealVar('N1%i'%x, 'N1%i'%x, 0.8, 0, 1)
    # N2 = root.RooRealVar('N2%i'%x, 'N2%i'%x, 0.1*1, 0, 0.2)
    params.extend([N1,N2,N3])
    models.append(build_3exp([N1,N2,N3], 'mjjBinModel%i'%x))
    simult.addPdf(models[-1], 'mjjBin%i'%x)



# fit
k = root.RooRealVar('k', 'k', -.1, -1, 0)
# q = root.RooRealVar('q', 'q', -.1, -10, 0)
exp1 = root.RooExponential('dummy1', 'dummy1', met, k)
# exp2 = root.RooExponential('dummy2', 'dummy2', met, q)
# norm = root.RooRealVar('dummy3', 'dummy3', 0.5*integral, 0, integral)
# norm2 = root.RooRealVar('dummy4', 'dummy4', 0.5*integral, 0, integral)
# dummy_model = root.RooAddPdf('add','add',root.RooArgList(exp1,exp2), root.RooArgList(norm,norm2))
# models[0].Print('v')
# models[0] = exp1
# models[0].fitTo(dh_data[0], root.RooFit.SumW2Error(True))
# models[0].fitTo(dh_data[0], root.RooFit.SumW2Error(True), root.RooFit.Extended(), root.RooFit.Strategy(2))
# models[0].fitTo(dh_data[0], root.RooFit.SumW2Error(True), root.RooFit.Extended())


# dsample = root.RooCategory('dsample','')
# dsample.defineType('pass')
# dsample.defineType('fail')
# dsimult = root.RooSimultaneous('dsimult','dsimult',dsample)
# dsimult.addPdf(models[0],'pass')
# dsimult.addPdf(models[0],'fail')
# ddh = root.RooDataHist('dd','dd',root.RooArgList(met),root.RooFit.Index(dsample),root.RooFit.Import('pass',dh_data[0]),root.RooFit.Import('fail',dh_data[0]))
# dsimult.fitTo(ddh, root.RooFit.SumW2Error(True))

fitresult = simult.fitTo(data, 
                         root.RooFit.Extended(),
                         root.RooFit.Strategy(2),
                         root.RooFit.NumCPU(4),
                         root.RooFit.SumW2Error(True))

# plot 
plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.Logy(True)
plot.InitLegend()
root.gStyle.SetOptStat(0)

# parameters = {x:[] for x in xrange(NPARAM)}
# paramerrors = {x:[] for x in xrange(NPARAM)}

for x in xrange(NBINSX):
    plot.Reset(True)
    plot.AddCMSLabel()
    plot.AddSqrtSLabel()
    xaxis = h2.GetXaxis()
    lo = xaxis.GetBinLowEdge(x+1)
    hi = lo + xaxis.GetBinWidth(x+1)
    plot.AddPlotLabel('%.0f < m_{jj} < %.0f [GeV]'%(lo, hi), .18, .8, False, 42, .04)
    plot.AddHistogram(h_data[x], 'Monte Carlo', root.kData)

    # print k1.getVal()
    tfit = models[x].asTF(root.RooArgList(met))
    h = h_data[x]
    scale = sum([h.GetBinContent(y)*h.GetBinWidth(y) for y in xrange(1,NBINSY+1)]) / tfit.Integral(50,1000)
    print scale
    print  sum([h.GetBinContent(y)*h.GetBinWidth(y) for y in xrange(1,NBINSY+1)])
    # def f(x):
    #     return scale * tfit.Eval(float(x[0]))
    tfit2 = root.TF1('scaled_%i'%x, lambda x : scale * tfit.Eval(float(x[0])), 50, 1000)
    tfit2.SetLineColor(2)
    print tfit2.Integral(50., 1000.), h_data[x].Integral()
    plot.AddAdditional(tfit2, 'l', 'Fit')
    # plot.AddAdditional(double_exp, 'l', 'Fit')
    # for x in xrange(NPARAM):
    #     parameters[x].append(double_exp.GetParameter(x))
    #     paramerrors[x].append(double_exp.GetParError(x))

    plot.Draw(args.outdir, 'fit_simult_%i'%x)

# plot.Logy(False)
# for x in xrange(NPARAM):
#     plot.Reset(True)
#     plot.AddCMSLabel()
#     plot.AddSqrtSLabel()
#     h_p = h_params[x]
#     for iMjj in xrange(0, h.GetNbinsX()):
#         h_p.SetBinContent(iMjj+1, parameters[x][iMjj])
#         h_p.SetBinError(iMjj+1, paramerrors[x][iMjj])
#     plot.AddHistogram(h_p, 'Fit parameter '+double_exp.GetParName(x), root.kSignal1, root.kSignal1, 'elp')
#     plot.Draw(args.outdir,'fit_param%i'%x)
