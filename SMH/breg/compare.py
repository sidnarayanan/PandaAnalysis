#!/usr/bin/env python

from sys import argv

v0 = argv[1]
v1 = argv[2]
binning = '>>h%i(50,'+','.join(argv[3:])+')' if len(argv)>3 else '>>h%i'

argv = []

import ROOT as root
root.gROOT.SetBatch(1)
from PandaCore.Utils.load import Load
Load('CanvasDrawer')

plot = root.CanvasDrawer()
plot.SetTDRStyle("vbf")
root.gStyle.SetOptStat(1)
root.gStyle.SetCanvasDefH(800)
root.gStyle.SetCanvasDefW(1200)

f0 = root.TFile.Open('/data/t3home000/snarayan/ttbar.root')
f1 = root.TFile.Open('/mnt/hadoop/scratch/snarayan/panda/breg_010_0test/batch/TTToHadronic_CP5_0_0.root')
#f1 = root.TFile.Open('/tmp/snarayan/TTTo2L2Nu_CP5_6_0.root')

t0 = f0.Get('jets')
t1 = f1.Get('events')

c0 = root.TCanvas()
plot.SetCanvas(c0)
c0.Divide(2)

v1_ = 'max(0,%s)'%(v1)
t1.Draw(v1_+binning%1,'jotFlav==5&&fabs(jotEta)<2.5','goff')

v0_ = 'max(0,%s)'%(v0)
t0.Draw(v0_+binning%0,'','goff')

c0.cd(1)
root.h0.SetName('theirs')
root.h0.SetLineWidth(2)
root.h0.SetNormFactor()
root.h0.GetXaxis().SetTitle(v0)
root.h0.SetTitle('')
root.h0.Draw()

c0.cd(2)
root.h1.SetName('ours')
root.h1.SetLineWidth(2)
root.h1.SetNormFactor()
root.h1.GetXaxis().SetTitle(v1)
root.h1.SetTitle('')
root.h1.Draw()

outpath = '/home/snarayan/public_html/figs/smh/breg/compare/'
v0 = v0.replace('*','').replace('/','')
v1 = v1.replace('*','').replace('/','')
plot.Draw(outpath,'%s_%s'%(v0,v1))

c0.cd(1).SetLogy()
c0.cd(2).SetLogy()
plot.Draw(outpath,'%s_%s_logy'%(v0,v1))

'''
import readline 
import code
variables = globals().copy()
variables.update(locals())
shell = code.InteractiveConsole(variables)
shell.interact()
'''
