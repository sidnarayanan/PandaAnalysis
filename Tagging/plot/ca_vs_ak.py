#!/usr/bin/env python

from PandaCore.Tools.script import * 
from PandaCore.Tools.root_interface import Selector 
import numpy as np

Load('HistogramDrawer')
args = parse('--outdir')

base = {}
base['ca'] = flatdir.replace('_ak_','_ca_')
base['ak'] = base['ca'].replace('_ca_','_ak_')

algos = ['ca','ak']
procs = ['QCD', 'ZpTT_lo']
cuts = {}
cuts['QCD'] = 'gen_pt>200 && fabs(gen_eta)<2.5'
cuts['ZpTT_lo'] = tAND(cuts['QCD'], 'fabs(gen_pdgid)==6 && gen_size<2.25')
weight = 'ptweight'
branches = ['gen_pt','clf_IsMatched', 'clf_MSD', 'clf_M', 'clf_Pt', 'clf_Tau32SD', 'clf_Tau32', weight]

plot = root.HistogramDrawer()
#plot.Stack(True)
#plot.Ratio(True)
#plot.FixRatio(0.4)
plot.SetTDRStyle()
plot.InitLegend(.2, .68, .65, .9, 2)
plot.DrawMCErrors(True)

sel = {}
for algo in algos:
    sel[algo] = {}
    for proc in procs:
        sel[algo][proc] = Selector()
        sel[algo][proc].read_files([base[algo] + '/' + proc + '.root'],
                                   branches=branches, cut=cuts[proc])



def draw(field, hbase, xlabel):
    h = {algo:{proc:sel[algo][proc].draw(fields=field, weight=weight, hbase=hbase)
               for proc in procs}
         for algo in algos}

    i = 1
    for ip,proc in enumerate(procs):
        for ia,algo in enumerate(algos):
            h_ = h[algo][proc]
            h_.SetLineStyle(ia+1)
            h_.Scale(1./h_.Integral())
            h_.GetXaxis().SetTitle(xlabel)
            plot.AddHistogram(h_, '%s %s'%(proc.replace('ZpTT_lo','Top'), algo.upper()), i, ip+1)
            i += 1
        
    plot.Draw(args.outdir+'/', field)
    plot.Reset()


draw('clf_Tau32', root.TH1D('','',50,0,1), '#tau_{32}')
draw('clf_Tau32SD', root.TH1D('','',50,0,1), '#tau_{32}^{SD}')
draw('clf_MSD', root.TH1D('','',50,0,400), 'Jet m_{SD} [GeV]')
draw('clf_M', root.TH1D('','',50,0,400), 'Jet m [GeV]')
draw('clf_Pt', root.TH1D('','',50,200,1000), 'Jet p_{T} [GeV]')
draw('gen_pt', root.TH1D('','',50,200,1000), 'Parton p_{T} [GeV]')
