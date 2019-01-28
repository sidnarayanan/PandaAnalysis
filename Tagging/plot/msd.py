#!/usr/bin/env python

from PandaCore.Tools.script import * 
from PandaCore.Tools.root_interface import Selector 
import numpy as np

Load('HistogramDrawer')
args = parse('--outdir')

procs = ['QCD', 'ZpTT_lo']
cuts = {}
cuts['QCD'] = 'gen_pt>200 && fabs(gen_eta)<2.5'
cuts['ZpTT_lo'] = tAND(cuts['QCD'], 'fabs(gen_pdgid)==6 && gen_size<1.44')
weight = 'ptweight'
branches = ['gen_pt','clf_MSD', 'clf_M', 'clf_Tau32', 'clf_Tau32SD', weight]
pt_bins = [200, 300, 500, 700, 1000]
#colors = [root.kBlack, root.kGray+3, root.kGray+2, root.kGray+1, root.kGray]
colors = [root.TColor.GetColor(r,0,b) for r,b in zip(np.linspace(1,0,len(pt_bins)),
                                                     np.linspace(0,1,len(pt_bins)))]
print colors 
#pt_bins.insert(0, 'Inclusive')

plot = root.HistogramDrawer()
#plot.Stack(True)
#plot.Ratio(True)
#plot.FixRatio(0.4)
plot.SetTDRStyle()
plot.InitLegend(.4, .7, .9, .9, 2)
#plot.InitLegend(.7, .3, .9, .9, 1)
plot.DrawMCErrors(True)
#plot.Logy()
#plot.SetAbsMin(5e-3)

sel = {}
for proc in procs:
    sel[proc] = Selector()
    sel[proc].read_files([flatdir + '/' + proc + '.root'],
                          branches=branches, cut=cuts[proc])



hbase = root.TH1D('','',50,50,450)

def draw(field, xlabel, proc):
    s = sel[proc]
    inclusive = None
    arg = []
    for i,b in enumerate(pt_bins):
        label = 'p_{T}>%i'%b
        mask = s['gen_pt']>b
        h = s.draw(fields=field, weight=weight, hbase=hbase, mask=mask)
        if inclusive is None:
            inclusive = h.Integral()
        h.Scale(1./h.Integral())
        h.SetLineColorAlpha(1, 1./(i+1))
#        h.SetMaximum(2e-1)
        h.GetXaxis().SetTitle(xlabel)
        arg.append((h, label, i+4, colors[i]))

    arg.reverse()
    for a in arg:
        plot.AddHistogram(*a)
        
    plot.Draw(args.outdir+'/', 'norm_'+field+'_'+proc)
    plot.Reset()


#draw('clf_M', 'Jet m [GeV]', 'QCD')
#draw('clf_M', 'Jet m [GeV]', 'ZpTT_lo')
#draw('clf_MSD', 'Jet m_{SD} [GeV]', 'QCD')
#draw('clf_MSD', 'Jet m_{SD} [GeV]', 'ZpTT_lo')
#
#hbase = root.TH1D('','',50,0,1)
#draw('clf_Tau32', '#tau_{32}', 'QCD')
#draw('clf_Tau32', '#tau_{32}', 'ZpTT_lo')
#draw('clf_Tau32SD', '#tau_{32}^{SD}', 'QCD')
#draw('clf_Tau32SD', '#tau_{32}^{SD}', 'ZpTT_lo')

hbase = root.TH1D('','',50,50,450)
def compare(proc):
    s = sel[proc]
    for field, label, color in [('clf_M', 'Mass', colors[-1]),
                                ('clf_MSD', 'SD Mass', colors[0])]:
        h = s.draw(fields=field, weight=weight, hbase=hbase)
        h.Scale(1./h.Integral())
        h.GetXaxis().SetTitle('Jet mass [GeV]')
        plot.AddHistogram(h, label, 5, color)
    plot.Draw(args.outdir+'/', 'compare_'+proc)
    plot.Reset()


compare('ZpTT_lo')
