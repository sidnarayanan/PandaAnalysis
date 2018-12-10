#!/usr/bin/env python

from PandaCore.Tools.script import * 
import numpy as np
from PandaCore.Tools.root_interface import Selector

Load('CanvasDrawer')
c = root.CanvasDrawer(800, 600)
c.SetTDRStyle()
c.InitLegend(.55,.7,.93,.93)
c.GetCanvas().SetRightMargin(0.05)
c.GetCanvas().SetTopMargin(0.05)
c.GetCanvas().SetLeftMargin(0.15)
c.GetCanvas().SetBottomMargin(0.15)
c.cd()

args = parse('--outdir', '--proc', '--field')

f = {}
f['pf'] = root.TFile.Open(flatdir + '/%s.root'%args.proc)
f['puppi'] = root.TFile.Open(flatdir.replace('nopuppi','ca') + '/%s.root'%args.proc)
t = {k:v.Get('events') for k,v in f.iteritems()}
sel = {k:Selector() for k in f}

is_top = args.proc == 'ZpTT_lo'
cut = 'gen_pt>200 && fabs(gen_eta)<2.5'
if is_top:
    cut = tAND(cut, 'gen_size<1.44 && fabs(gen_pdgid)==6') 

for k,s in sel.iteritems():
    s.read_tree(t[k], cut=cut, 
                branches=['npv', args.field, 'ptweight'])

nbins = 8
if args.field == 'clf_MSD':
    hbase = root.TH2D('','',nbins,1,40,50,0,400)
    hbase.GetYaxis().SetTitle('m_{SD} [GeV]')
else:
    hbase = root.TH2D('','',nbins,1,40,50,0.2,1.2)
    hbase.GetYaxis().SetTitle('#tau_{32}^{SD}')
hbase.GetXaxis().SetTitle('N_{PV}')
hbase.GetXaxis().SetTitleOffset(1); hbase.GetXaxis().SetTitleSize(0.06)
hbase.GetYaxis().SetTitleOffset(1); hbase.GetYaxis().SetTitleSize(0.06)

def draw(key, color):
    s = sel[key]
    h = s.draw(fields=['npv', args.field], hbase=hbase, weight='ptweight')
    h.SetMarkerColor(color)
    h.SetMarkerSize(1)
    h.SetFillColorAlpha(color, 0.1)
    h.SetLineColor(color)
    h.SetMarkerStyle(4)
    return h

h = {
        'pf' : draw('pf', root.kBlue),
        'puppi' : draw('puppi', root.kRed),
    }


#style = 'VIOLINX(03000130)'
style = 'CANDLEX(10000131)'
#style = 'violin'
h['pf'].Draw(style)
c.GetLegend().AddEntry(h['pf'], 'All particles', 'lfp')
h['puppi'].Draw(style+'SAME')
c.GetLegend().AddEntry(h['puppi'], 'PUPPI', 'lfp')
c.GetCanvas().Update()

c.Draw(args.outdir, 'npv_'+args.field+'_'+args.proc)

