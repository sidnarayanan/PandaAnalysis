#!/usr/bin/env python

from PandaCore.Tools.script import * 
import numpy as np
from PandaCore.Tools.root_interface import Selector 

baseDir = getenv('PANDA_FLATDIR')+'/' 
args = parse('--outdir', ('--cut', {'default':'1==1'}))

Load('PandaCoreDrawers')
plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.InitLegend()
# plot.Logy(True)

s = Selector()
f = root.TFile.Open(baseDir + '/GJets.root')
t = f.Get('events')
s.read_tree(t, 
            branches=['trueGenBosonPt', 
                      'normalizedWeight',
                      'normalizedWeight*sf_qcdV', 
                      'normalizedWeight*sf_qcdV*sf_ewkV'], 
            cut='trueGenBosonPt>160')

ptbins = (10,160,1600) 
hlo = s.draw('trueGenBosonPt', weight='normalizedWeight', fbins=ptbins)
hqcd = s.draw('trueGenBosonPt', weight='normalizedWeight*sf_qcdV', fbins=ptbins)
hewk = s.draw('trueGenBosonPt', weight='normalizedWeight*sf_qcdV*sf_ewkV', fbins=ptbins)

for h in [hewk, hqcd, hlo]:
  h.Divide(hlo)
for h in [hewk, hqcd, hlo]:
  h.Smooth()

hlo.GetXaxis().SetTitle('p_{T}^{Z} [GeV]')
hlo.GetYaxis().SetTitle('k')

plot.AddHistogram(hlo, '')
plot.AddHistogram(hqcd, 'k_{QCD}')
plot.AddHistogram(hewk, 'k_{QCD} #times k_{EWK}')

plot.Draw(args.outdir+'/', 'acorr_ptv')
