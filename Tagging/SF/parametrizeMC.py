#!/usr/bin/env python

from sys import argv,exit
argv=[]

import ROOT as root
from PandaCore.Utils.load import *
import numpy as np
from scipy.optimize import minimize, curve_fit

np.set_printoptions(linewidth='200')
Load('HistogramDrawer') 

basedir = '/home/snarayan/public_html/figs/toptagging/v9/sf/'
files = {
    'pass' : root.TFile.Open(basedir + '/pass_tight_hists.root')
    }

def gaus(x, mu, sigma, N):
  return N * np.exp(-np.power(x - mu, 2) / (sigma ** 2)) / np.sqrt(6.28 * sigma**2)

def gaus2(x, *args):
  return  gaus(x, *args[0:3]) + gaus(x, *args[3:6])

def gaus_lin(x, a, b, mu, sigma, N):
  return N * (a*x + b + np.exp(-np.power(x - mu, 2) / (sigma ** 2)) / np.sqrt(6.28 * sigma**2))

def exp(x, alpha, N):
  return N * np.exp(-alpha * x) 

plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.InitLegend()
root.gStyle.SetOptStat(0)

def fitter(hist, fn, init_params, name):
  # convert hist to array 
  nbins = hist.GetNbinsX()
  bins = range(1, nbins+1)
  bins = range(1, 30)
  x = np.array([hist.GetBinCenter(i) for i in bins])
  y = np.array([hist.GetBinContent(i) for i in bins])
  s = np.array([hist.GetBinError(i) for i in bins])

  fit_params, cov = curve_fit(fn, xdata=x, ydata=y, sigma=s, 
                              p0=init_params, absolute_sigma=True)

  x_plot = np.linspace(x[0], x[-1], 200)

  # sample some set of fit parameter points 
  sample_params = np.random.multivariate_normal(fit_params, cov, 100)
  zeros = np.zeros(x_plot.shape)
  yhat = np.zeros(x_plot.shape)
  yhat_hi = np.zeros(x_plot.shape)
  yhat_lo = np.zeros(x_plot.shape)
  for i, x_ in enumerate(x_plot):
    yhat[i] = fn(x_, *fit_params)
    yhat_dist = sorted([fn(x_, *p) for p in sample_params])
    yhat_lo[i] = yhat[i] - yhat_dist[16] 
    yhat_hi[i] = yhat_dist[84] - yhat[i]


  # convert fit result to TGraph
  g = root.TGraphAsymmErrors(len(x_plot), x_plot, yhat,
                             zeros, zeros,
                             yhat_lo, yhat_hi)
  g.SetLineColor(root.kRed)
  g.SetFillColorAlpha(root.kRed, 0.5)
  g.SetLineWidth(1)
  plot.Reset(False)
  plot.AddHistogram(hist, 'Template', root.kData)
  plot.AddAdditional(g, 'l3', 'Fit')
  plot.Draw(basedir, 'param_'+name)

  return fit_params, cov


p,c = fitter(files['pass'].Get('h_fjMSD_3Minusprong'), gaus2, ([175, 20, 5000] * 2), 'pass_3prong')

p,c = fitter(files['pass'].Get('h_fjMSD_1Minusprong'), exp, [10, 10000], 'pass_1prong')
p,c = fitter(files['pass'].Get('h_fjMSD_2Minusprong'), exp, [10, 10000], 'pass_2prong')
