#!/usr/bin/env python

from PandaCore.Tools.script import * 
from PandaCore.Drawers.plot_utility import *
import PandaAnalysis.Monotop.CombinedBVetoSelection as sel
import json 

args = parse('--outdir', ('--ecf', STORE_TRUE))

lumi = 35900.
cut = 'pt>250 && mSD>110 && mSD<210 && fabs(eta)<2.5'

plot = PlotUtility()
plot.n_threads = 10
plot.SetTDRStyle()
plot.InitLegend()
plot.cut = cut
plot.SetLumi(lumi/1000)
plot.SetNormFactor(True)
plot.AddSqrtSLabel()
plot.do_overflow = True
plot.do_underflow = True

plot.mc_weight = 'normalizedWeight'

lqg = Process('LQG jets', root.kWjets, root.kBlue+1)
lqg.additional_weight = 'ptweight'

top =  Process('Top jets',root.kTTbar, root.kRed+2)
top.additional_cut = 'matched==1 && gensize<1.2'
top.additional_weight = 'ptweight'

flatdir = '/home/snarayan/home000/store/scramjet/v7/'
lqg.add_file(flatdir + 'QCD_evt8.root')
top.add_file(flatdir + 'ZpTT.root')

for p in [lqg, top]:
    p.tree_name = 'puppiCA15'
    plot.add_process(p)

with open('ranges.json','r') as fp:
    payload = json.load(fp)
    ratios = payload['ratios']
    ecfs = payload['ecfs']

def p2s(label):
    return label.replace('fjECFN','ecfN').replace('[0]','')

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
        plot.add_distribution(FDistribution(p2s(e['formula']), e['lo'], e['hi'], 20, e['label'], 'a.u.',
                                            filename=filename))

    for i,r in enumerate(ratios):
        if r['hi'] <= r['lo']:
            continue
        f = eval(p2s(r['formula']))
        filename = 'ratio_'
        filename += convert(r['label'].split('/')[0])
        filename += convert(r['label'].split('/')[1].split('^')[0])
        plot.add_distribution(FDistribution(f,
                                            r['lo'], r['hi'], 20,
                                            r['label'], 'a.u.',
                                            filename=filename))

for a in [('tau32', 0, 1.2, 20, '#tau_{32}', 'a.u.'),
          ('tau32SD', 0, 1.2, 20, '#tau_{32}^{SD}', 'a.u.'),
          ('top_ecfv8_bdt', -1.2, 1.2, 20, 'Combined BDT', 'a.u.'),
          ('top_allv2_bdt', -1.2, 1.2, 20, '50 ECF BDT', 'a.u.'),
          ('top_ecfv6_bdt', -1.2, 1.2, 20, '11 ECF BDT', 'a.u.'),
          ('top_taufrec_bdt', -1.2, 1.2, 20, '#tau_{32}^{SD}+f_{rec} BDT', 'a.u.'),
          ('pt', 250,1000,20,'CA15 p_{T} GeV','a.u.'), 
          ]:
    plot.add_distribution(FDistribution(*a, filename=a[0].replace('[0]','')))

plot.draw_all(args.outdir+'/')
