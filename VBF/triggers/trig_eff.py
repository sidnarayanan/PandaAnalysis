#!/usr/bin/env python

from PandaCore.Tools.script import * 
from PandaCore.Tools.root_interface import Selector
from numpy import arange
from math import sqrt
import PandaAnalysis.VBF.PandaSelection as sel 

basedir = getenv('PANDA_FLATDIR') + ''
infile = basedir+'/SingleMuon.root'

args = parse('--outdir', ('--nmu', {'type':int, 'default':1}))
figsdir = args.outdir

Load('GraphAsymmErrDrawer')

lumi=36000
if args.nmu == 1:
  base_cut = sel.cuts['singlemuon']
  metnomu = 'pfUWmag'
else:
  base_cut = sel.cuts['dimuon']
  metnomu = 'pfUZmag'
base_cut = tAND(base_cut, '(trigger&8)!=0')
# base_cut = remove_cut(base_cut, 'pfUWmag')
# base_cut = tAND(base_cut, 'pfUWmag>180')
base_cut = base_cut.replace('pfUWmag>250', 'pfUWmag>200')

plot = root.GraphAsymmErrDrawer()
plot.SetTDRStyle()
plot.SetLumi(lumi/1000.)
plot.DrawEmpty(True)
#plot.SetAutoRange(False)
root.gStyle.SetOptStat(0)


s = Selector()
f = root.TFile.Open(infile); t = f.Get('events')
branches = [metnomu, 'jot12Mass', 'jetPt[0]', 'jetPt[0]+jetPt[1]', 'fabs(jotEta[0])', 'trigger&1', 'fabs(jotEta[1])', 'barrelHT', 'barrelHTMiss']
s.read_tree(t, branches = branches, cut = base_cut)

triggered = ( s['trigger&1'] != 0 )

j1 = (s['fabs(jotEta[0])'] < 3)
j2 = (s['fabs(jotEta[1])'] < 3)
bb = (j1 & j2).flatten()
bf = (j1 ^ j2).flatten()
ff = ~(j1 | j2).flatten()
bins = [200, 225, 250, 275, 350, 400, 450, 500, 550, 600, 650,  750]
ratios_to_plot = {
  metnomu : (bins,'U [GeV]'),
  'jot12Mass' : ([200, 400, 600, 800, 1000, 1250, 1500, 1750, 2000, 2500, 3500, 5000],'m_{jj} [GeV]'),
  'jetPt[0]' : ([0, 40, 80, 120, 160, 200, 300, 500],'Central jet 1 p_{T} [GeV]'),
  'jetPt[0]+jetPt[1]' : ([0, 40, 80, 120, 160, 200, 300, 500],'Central jet 1 p_{T} + jet 2 p_{T} [GeV]'),
#   'barrelJet1Pt' : ([0, 40, 80, 120, 160, 200, 300, 500],'Barrel jet 1 p_{T} [GeV]'),
#   'barrelJet12Pt' : ([0, 40, 80, 120, 160, 200, 300, 500],'Barrel jet 1 + jet 2 p_{T} [GeV]'),
  'barrelHT' : ([0, 80, 120, 160, 200, 240, 280, 320, 360, 400, 450, 500, 600, 1000],'Barrel H_{T} [GeV]'),
  'barrelHTMiss' : (bins,'Barrel H_{T}^{miss} [GeV]'),
}

fout = root.TFile(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/vbf16/trig/param_nmu%i.root'%args.nmu, 'RECREATE')

for k,v in ratios_to_plot.iteritems():
    plot.Clear()
    plot.InitLegend(0.4, 0.8, 0.95, 0.9, 3)
#    plot.AddCMSLabel()
    plot.AddLumiLabel()

    h_inc = s.draw(k, weight=None, vbins=v[0])
    h_trig = s.draw(k, weight=None, mask=triggered, vbins=v[0])

    hh_ratio = h_trig.Clone(); hh_ratio.Divide(h_inc)
    h_ratio = root.TGraphAsymmErrors()
    h_ratio.BayesDivide(h_trig, h_inc)
#    for ib in xrange(1, h_ratio.GetNbinsX()+1):
#      h_ratio.SetBinContent(ib, h_trig.GetBinContent(ib) / h_inc.GetBinContent(ib))
#      h_ratio.SetBinError(ib, h_trig.GetBinError(ib) / h_inc.GetBinContent(ib))
#    h_ratio.Divide(h_inc)
    h_ratio.GetXaxis().SetTitle(v[1])
    h_ratio.GetYaxis().SetTitle('Efficiency')
    h_ratio.GetYaxis().SetTitleOffset(1.4)
    if hh_ratio.GetMinimum() < 0.5:
      h_ratio.SetMaximum(1.4)
      h_ratio.SetMinimum(0)
    else:
      height = hh_ratio.GetMaximum() - hh_ratio.GetMinimum()
      h_ratio.SetMaximum(1 + height * 0.4)
      h_ratio.SetMinimum(hh_ratio.GetMinimum() - (height * 0.4))

#     fout.WriteTObject(hh_ratio, 'h_%s'%k)
#     fout.WriteTObject(h_ratio, 'g_%s'%k)

    plot.AddGraph(h_ratio, 'All', 1, 1, 'ep')

    h_ratio_masked = {}
    for i, (mask, label) in enumerate([(bb, 'BB'), (bf, 'BF')]):
        h_inc = s.draw(k, weight=None, mask=mask, vbins=v[0])
        h_trig = s.draw(k, weight=None, mask=(mask & triggered), vbins=v[0])
        h_ratio_masked[label] = h_ratio.Clone() 
        h_ratio_masked[label].BayesDivide(h_trig, h_inc)
        hh_ratio = h_trig.Clone(); hh_ratio.Divide(h_inc)
        fout.WriteTObject(hh_ratio, 'h_%s_%s'%(k,label))
        fout.WriteTObject(h_ratio_masked[label], 'g_%s_%s'%(k,label))

        plot.AddGraph(h_ratio_masked[label], label, i+2, 1, 'ep')

    plot.Draw(args.outdir,'trigeff_nmu%i'%args.nmu+k)
