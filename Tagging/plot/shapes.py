#!/usr/bin/env python

from PandaCore.Tools.script import * 
from PandaCore.Drawers.plot_utility import *
import PandaAnalysis.Monotop.CombinedBVetoSelection as sel
import json 

args = parse('--outdir')

lumi = 35900.
region = 'singlemuonw' 
cut = sel.cuts[region]
weight = sel.weights[region]%lumi

plot = PlotUtility()
plot.Stack(True)
plot.Ratio(True)
plot.FixRatio(0.4)
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.cut = cut
plot.SetLumi(lumi/1000)
plot.SetNormFactor(True)
plot.AddLumiLabel(True)
plot.do_overflow = True
plot.do_underflow = True

weight = sel.weights[region]%lumi
plot.mc_weight = weight

wjets         = Process('W+jets',root.kWjets)
diboson        = Process('Diboson',root.kDiboson)
qcd             = Process("QCD",root.kQCD)
top            = Process('Top',root.kTTbar)
data            = Process("Data",root.kData)
processes = [qcd,diboson,top,wjets,data]

wjets.add_file(flatdir+'WJets.root')
top.add_file(flatdir+'TTbar.root')
top.add_file(flatdir+'SingleTop.root')
diboson.add_file(flatdir+'Diboson.root')
qcd.add_file(flatdir+'QCD.root')
data.add_file(flatdir+'MET.root')
data.additional_cut = sel.triggers['met']

for p in processes:
    plot.add_process(p)

with open('ranges.json','r') as fp:
    payload = json.load(fp)
    ratios = payload['ratios']
    ecfs = payload['ecfs']

for i,e in enumerate(ecfs):
    if e['hi'] <= e['lo']:
        continue 
    plot.add_distribution(FDistribution(e['formula'], e['lo'], e['hi'], 20, e['label'], 'a.u.',
                                        filename='ecf_%i'%i))

for i,r in enumerate(ratios):
    if r['hi'] <= r['lo']:
        continue
    f = eval(r['formula'])
    plot.add_distribution(FDistribution(f,
                                        r['lo'], r['hi'], 20,
                                        r['label'], 'a.u.',
                                        filename='ratio_%i'%i))

plot.draw_all(args.outdir+'/'+region+'_')
