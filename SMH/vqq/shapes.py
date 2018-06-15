#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
dataDir = baseDir

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cat',metavar='cat',type=str,default='resolved')
args = parser.parse_args()

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.SMH.VqqSel as sel
from PandaCore.Drawers.plot_utility import *

cut = sel.cuts['signal']

### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.Stack(False)
plot.SetTDRStyle()
plot.InitLegend()
plot.AddCMSLabel()
plot.SetNormFactor(True)
plot.AddSqrtSLabel()
plot.cut = cut
plot.do_overflow = True
plot.do_underflow = True

weight = sel.weights['signal']%1
plot.mc_weight = weight

### DEFINE PROCESSES ###
vh           = Process('V(#rightarrowqq)H(#rightarrowbb)',root.kExtra4)
qcd           = Process('QCD', root.kQCD)

### ASSIGN FILES TO PROCESSES ###
vh.add_file(baseDir+'VqqHbb.root')
qcd.add_file(baseDir+'QCD.root')

processes = [qcd, vh]

for p in processes:
    plot.add_process(p)

def add(*args, **kwargs):
    plot.add_distribution(FDistribution(*args, **kwargs))

add('fjMSD[fjHiggsIdx]', 0, 300, 30, 'm_{H} [Gev]', 'a.u.', filename='fjHMSD')
add('fjMSD[fjVIdx]', 0, 300, 30, 'm_{V} [Gev]', 'a.u.', filename='fjVMSD')
add('fjMSD_corr[fjHiggsIdx]', 0, 300, 30, 'Corr m_{H} [Gev]', 'a.u.', filename='fjHMSD_corr')
add('fjMSD_corr[fjVIdx]', 0, 300, 30, 'Corr m_{V} [Gev]', 'a.u.', filename='fjVMSD_corr')
add('fjEta[fjHiggsIdx]', -2.5, 2.5, 30, '#eta_{H}', 'a.u.', filename='fjHEta')
add('fjEta[fjVIdx]', -2.5, 2.5, 30, '#eta_{V}', 'a.u.', filename='fjVEta')
add('fjPt[fjHiggsIdx]', 250, 1000, 30, 'p_{T,H} [Gev]', 'a.u.', filename='fjHPt')
add('fjPt[fjVIdx]', 250, 1000, 30, 'p_{T,V} [Gev]', 'a.u.', filename='fjVPt')
add('fabs(SignedDeltaPhi(fjPhi[0], fjPhi[1]))', 0, 3.1416, 30, '|#Delta_{#phi}(H,V)|', 'a.u.', filename='fjDPhi')
add('(fjPt[fjHiggsIdx] - fjPt[fjVIdx]) / fjPt[fjHiggsIdx]', -1, 1,  30, '(p_{T,H} - p_{T,V}) / p_{T,H}', 'a.u.', filename='fjPtBalance')
add(sel.n2('fjHiggsIdx'), 0, 0.5, 30, 'N2(H)', 'a.u.', filename='fjHN2')
add(sel.n2('fjVIdx'), 0, 0.5, 30, 'N2(V)', 'a.u.', filename='fjVN2')
add('fjDoubleCSV[fjHiggsIdx]', 0, 1, 30, 'DoubleCSV(H)', 'a.u.', filename='fjHDoubleCSV')
add('fjDoubleCSV[fjVIdx]', 0, 1, 30, 'DoubleCSV(V)', 'a.u.', filename='fjVDoubleCSV')
add('(fjDoubleCSV[fjHiggsIdx] - fjDoubleCSV[fjVIdx])', 0, 2,  30, 'DoubleCSV(H) - DoubleCSV(V)', 'a.u.', filename='fjDoubleCSVDifference')


### DRAW AND CATALOGUE ###
plot.draw_all(args.outdir+'/')
