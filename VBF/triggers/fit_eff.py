#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange
import numpy as np

basedir = getenv('PANDA_FLATDIR') + ''
infile = basedir+'/SingleMuon.root'


parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--infile',metavar='infile',type=str,default=None)
parser.add_argument('--variable',metavar='variable',type=str,default='barrelHTMiss')
args = parser.parse_args()

figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Tools.Load import Load
from PandaCore.Tools.Misc import *
from math import sqrt
Load('GraphAsymmErrDrawer')
root.gROOT.SetBatch()

lumi=36000
plot = root.GraphAsymmErrDrawer()
plot.SetTDRStyle()
plot.SetLumi(lumi/1000.)
plot.DrawEmpty(True)
plot.InitLegend(0.55,0.76,0.88,0.9)
#plot.SetAutoRange(False)
root.gStyle.SetOptStat(0)


fin = root.TFile(args.infile)
fout = root.TFile(args.infile.replace('param','fit'), 'RECREATE')
#ffit = root.TF1("sigmoidFit", "[2]*pow(1+[3]*TMath::Exp(-[0]*(x-[1])),-1./[4])", 10, 1000)
#ffit.SetParNames("k", "E_{c}","#varepsilon_{obs}","Q","#nu");
#ffit.SetParameter(0,.03);  
#ffit.SetParameter(1,160); 
#ffit.SetParameter(2,.95);
#ffit.SetParameter(3,1); 
#ffit.SetParameter(4,1);
#ffit.SetParLimits(0,0,100)
##ffit.SetParLimits(1,0,100)
#ffit.SetParLimits(2,0,1)
#ffit.SetParLimits(3,0,10)
#ffit.SetParLimits(4,0,100)
ffit = root.TF1("sigmoidFit", "[0]*TMath::Erf((x-[1])/[2]) + [3]*TMath::Erf((x-[4])/[5])", 0, 1000)
ffit.SetParameters(0.5, 100, 100, 0.5, 150, 150)
ffit.SetParLimits(0, 0, 1)
ffit.SetParLimits(3, 0, 1)


g = fin.Get('g_%s'%(args.variable))
result = g.Fit(ffit, 'S', '', 10, 1000)
result = g.Fit(ffit, 'S', '', 10, 1000)
prob = root.TMath.Prob(result.Chi2(), g.GetN() - 6)
print prob

fout.WriteTObject(g)
fout.WriteTObject(ffit, 'f_eff')
fout.Close()

plot.AddGraph(g, 'Efficiency', 2, 1, 'EP')
plot.AddAdditional(ffit, 'l', 'Fit, P(#chi^{2}|NDF) = %.2f'%prob)
outname = args.infile.split('/')[-1].replace('.root','')
plot.Draw(args.outdir,'fit_'+outname+'_'+args.variable)

