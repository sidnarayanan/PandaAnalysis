#!/usr/bin/env python

from PandaCore.Tools.script import * 
from PandaCore.Drawers.plot_utility import *
import json 

args = parse('--outdir', ('--masscut', STORE_TRUE))

plot = PlotUtility()
plot.SetTDRStyle()
plot.InitLegend(.65, .7, .93, .91)
plot.cut = Cut('filter_maxRecoil>175 && fjPt[0]>250 && fabs(fjEta[0])<2.5')
if args.masscut:
    plot.cut &= 'fabs(fjMSD[0]-160)<50'
    plot.AddPlotLabel('110 < m_{SD} < 210 GeV',.18,.85,False,42,.04)
plot.SetNormFactor(True)
plot.do_overflow = True
plot.do_underflow = True
plot.mc_weight = 'normalizedWeight'

qcd             = Process("LQG",root.kWjets, root.kBlue+1)
qcd.add_file(flatdir+'/WJets.root')

top            = Process('Top',root.kTTbar, root.kRed+2)
top.add_file(flatdir+'/ZpTT.root')
top.additional_cut = 'fjIsMatched[0]==1 && fjGenSize[0]<1.44'

for p in [qcd, top]:
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

'''
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
'''

for a in [('fjTau32[0]', 0, 1.2, 20, '#tau_{32}', 'a.u.'),
          ('fjTau32SD[0]', 0, 1.2, 20, '#tau_{32}^{SD}', 'a.u.'),
          ('fjHTTMass[0]', 0, 300, 20, 'HTT mass [GeV]', 'a.u.'),
          ('fjMSD[0]', 0, 300, 20, 'SD mass [GeV]', 'a.u.'),
          ('top_ecf_bdt', -1.2, 1.2, 20, 'Top ID BDT score', 'a.u.'),
          ('fjHTTFRec[0]', 0, 0.8, 20, 'HTT f_{rec}', 'a.u.')]:
    plot.add_distribution(FDistribution(*a, filename=a[0].replace('[0]','')))

name = '/mass_' if args.masscut else '/incl_'
plot.draw_all(args.outdir+name)
