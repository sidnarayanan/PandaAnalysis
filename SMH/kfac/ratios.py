#!/usr/bin/env python

from PandaCore.Tools.script import * 
basedir = getenv('PANDA_FLATDIR') + '/' 

args = parse('--outdir')
figsdir = args.outdir

from PandaCore.Tools.root_interface import Selector
Load('HistogramDrawer')

plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.AddSqrtSLabel()
plot.AddCMSLabel()
plot.InitLegend()

#cut = 'genJet0Pt>80 && genJet2Pt>30 && trueGenBosonPt>100 && lheHT>100'
cut = '1==1'
weight = 'normalizedWeight'
pt = '0.001*trueGenBosonPt'

s = {}
for a in ['lo_incl','nlo']:
    s[a] = Selector()
    f = root.TFile.Open((basedir + 'WJets_' + a +'.root'))
    t = f.Get('events')
    s[a].read_tree(t, branches = [pt, weight], cut = cut)
    print 'done with',a

print 'done reading'

def draw(num, den, hbase, label, color):
    x = '0.001*trueGenBosonPt'
    hnum = s[num].draw([x], weight=weight, hbase=hbase)
    hden= s[den].draw([x], weight=weight, hbase=hbase)
    hnum.Divide(hden)
    hnum.GetXaxis().SetTitle('V p_{T} [TeV]')
    hnum.GetYaxis().SetTitle('Ratio')
    plot.AddHistogram(hnum, label, 11, color, 'e hist') 
    

bins = array('f', [0.001 * x for x in [-100,0,30,60,90,120,150,180,210,240,280,320,360,400,450,500,550,600,700,800]])

hbase = root.TH1D('', '', len(bins)-1, bins)
draw('nlo', 'lo_incl', hbase, 'incl nlo/lo', root.kRed)

plot.Draw(figsdir, '/genbosonpt')
