#!/usr/bin/env python

from PandaCore.Tools.script import * 
from PandaCore.Drawers.plot_utility import *
import PandaAnalysis.Monotop.CombinedBVetoSelection as sel
import json 

args = parse('--outdir', '--region', ('--ecf', STORE_TRUE))

lumi = 35900.
region = args.region 
cut = sel.cuts[region]
weight = sel.weights[region]%lumi

plot = PlotUtility()
plot.Stack(True)
plot.Ratio(True)
plot.FixRatio(0.4)
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.cut = cut
plot.SetLumi(lumi/1000)
plot.SetNormFactor(True)
plot.AddLumiLabel(True)
plot.do_overflow = True
plot.do_underflow = True

weight = sel.weights[region]%lumi
plot.mc_weight = weight

zjets         = Process('Z+jets',root.kZjets)
wjets         = Process('W+jets',root.kWjets)
diboson        = Process('Diboson',root.kDiboson)
qcd             = Process("QCD",root.kQCD)
top            = Process('Top',root.kTTbar)
data            = Process("Data",root.kData)
if region == 'dimuon':
    processes = [diboson,top,zjets,data]
else:
    processes = [qcd,diboson,top,wjets,data]

zjets.add_file(flatdir+'ZJets.root')
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

def convert(label):
    label = label.strip('(')
    label = label.replace('e(','').replace(')','')
    label = label.split(',')
    return '%i%i%.2i'%(int(label[0]), int(label[1]), int(10*float(label[2])))

if args.ecf:
    for i,e in enumerate(ecfs):
        if e['hi'] <= e['lo']:
            continue 
        filename = 'ecf_' + convert(e['label'])
        plot.add_distribution(FDistribution(e['formula'], e['lo'], e['hi'], 20, e['label'], 'a.u.',
                                            filename=filename))

    for i,r in enumerate(ratios):
        if r['hi'] <= r['lo']:
            continue
        f = eval(r['formula'])
        filename = 'ratio_'
        filename += convert(r['label'].split('/')[0])
        filename += convert(r['label'].split('/')[1].split('^')[0])
        plot.add_distribution(FDistribution(f,
                                            r['lo'], r['hi'], 20,
                                            r['label'], 'a.u.',
                                            filename=filename))

for a in [('fjTau32[0]', 0, 1.2, 20, '#tau_{32}', 'a.u.'),
          ('fjTau32SD[0]', 0, 1.2, 20, '#tau_{32}^{SD}', 'a.u.'),
          ('fjHTTMass[0]', 0, 300, 20, 'HTT mass [GeV]', 'a.u.'),
          ('fjMSD[0]', 0, 300, 20, 'SD mass [GeV]', 'a.u.'),
          ('top_ecf_bdt', -1.2, 1.2, 20, 'Top ID BDT score', 'a.u.'),
          ('fjHTTFRec[0]', 0, 0.8, 20, 'HTT f_{rec}', 'a.u.')]:
    plot.add_distribution(FDistribution(*a, filename=a[0].replace('[0]','')))

plot.draw_all(args.outdir+'/'+region+'_')
