#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
args = parser.parse_args()
sname = argv[0]

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Utils.load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
from PandaCore.Drawers.plot_utility import * 

### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.Logy(True)
plot.Stack(True)
plot.DrawEmpty(True)
plot.Ratio(True)
plot.FixRatio(0.4)
plot.SetRatioLabel('#frac{NNLO-NLO}{NLO}')
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()

plot.mc_weight = 'normalizedWeight'
plot.cut = '1==1'

### DEFINE PROCESSES ###
nlo     = Process('NLO Powheg',root.kTTbar); 
processes = [nlo]

### ASSIGN FILES TO PROCESSES ###
nlo.add_file(baseDir+'TTbar.root')

for p in processes:
  plot.add_process(p)

ptbins = [0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,800,900,1000,1100,1200,1400,1600,1800,2000,2200,2600,3000]
nnlo_vals = [131.27,257.3,207.66,117.9,57.799,27.369,13.112,6.4831,3.3328,1.7778,0.98056,0.56038,0.32985,0.19885,0.1992,0.081692,0.036048,0.016542,0.0080277,0.0059894,0.0016044,0.00046105,0.00014635,4.6246e-05,2.071e-05,2.4198e-06] # mT/2
#nnlo_vals = [124.78,246.57,199.13,112.8,55.208,26.099,12.492,6.171,3.1758,1.695,0.93597,0.53512,0.31533,0.19045,0.19136,0.078658,0.034769,0.01601,0.0077749,0.0058412,0.0015839,0.0004631,0.00014697,4.6838e-05,2.1337e-05,2.5606e-06] #mT
h_nnlo = root.TH1D('hnnlo','',len(ptbins)-1,array('f',ptbins))
for ib in xrange(1,len(ptbins)):
    h_nnlo.SetBinContent(ib,nnlo_vals[ib-1])
    h_nnlo.SetBinError(ib,0.00001)
divide_bin_width(h_nnlo)
h_nnlo.GetXaxis().SetTitle('Average Top p_{T} [GeV]')
h_nnlo.GetYaxis().SetTitle('d#sigma/dp_{T}')

plot.add_distribution(VDistribution('(genTopPt+genAntiTopPt)/2',ptbins,'Average Top p_{T} [GeV]','d#sigma/dp_{T}',filename='aveTopPt'))

plot.canvas.AddHistogram(h_nnlo,'NNLO calc',root.kData)

plot.draw_all(args.outdir)
