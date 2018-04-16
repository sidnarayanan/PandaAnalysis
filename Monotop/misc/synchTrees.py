#!/usr/bin/env python

from sys import argv 
import argparse

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--sample',metavar='sample',type=str)
args = parser.parse_args()
argv = []

import ROOT as root
from PandaCore.Tools.Load import Load

Load('PandaCoreTools')
Load('PandaCoreDrawers')

es = root.EventSyncher()
es.relative = False

f1 = '/home/snarayan/home000/store/panda/v_8026_0_5_slim/%s.root'%args.sample
f2 = '/home/snarayan/home000/store/panda/v_8026_0_5_ak8/%s.root'%args.sample
es.preselection = 'fj1Pt>400'

formulae = [
		('fj1Pt',400,1200,'p_{T}',True),
		('fj1MaxCSV',0,1,'Max CA15 subjet CSV',False),
		]

for f in formulae:
	es.AddFormula(*f)

es.RunFiles(f1,f2,'events')

hists = es.PlayViolin(1)

plot = root.CanvasDrawer()
plot.SetTDRStyle()
c = root.TCanvas()

for iF in xrange(len(formulae)):
	hist = hists[iF]
	c.Clear()
	c.SetGrid()
	hist.SetFillColor(0);
	hist.SetLineColor(0);
	hist.Draw('VIOLIN')
	hist.GetYaxis().SetTitleOffset(1.75)
	if formulae[iF][4]:
		hist.GetYaxis().SetTitle('#LT(CA15 - AK8)/CA15#GT')
	else:
		hist.GetYaxis().SetTitle('#LTCA15 - AK8#GT')
	plot.AddCMSLabel(.16,.94)
	plot.SetLumi(35.8); plot.AddLumiLabel(True)
	plot.SetCanvas(c)
	plot.Draw(args.outdir+'/%s_'%args.sample,formulae[iF][0])

for iF in xrange(len(formulae)):
	hist = hists[iF]
	c.Clear()
	c.SetGrid()
	hist.SetFillColor(0);
	hist.SetLineColor(0);
	hist.Draw('VIOLIN')
	hist.GetYaxis().SetTitleOffset(1.75)
	if formulae[iF][4]:
		hist.GetYaxis().SetTitle('#LT(CA15 - AK8)/CA15#GT')
	else:
		hist.GetYaxis().SetTitle('#LTCA15 - AK8#GT')
	plot.AddCMSLabel(.16,.94)
	plot.SetLumi(35.8); plot.AddLumiLabel(True)
	plot.SetCanvas(c)
	plot.Draw(args.outdir+'/%s_'%args.sample,formulae[iF][0])
