#!/usr/bin/env python

from PandaCore.Tools.script import * 
import numpy as np
from PandaCore.Tools.root_interface import Selector

Load('CanvasDrawer')
c = root.CanvasDrawer(800, 600)
c.SetTDRStyle()
c.InitLegend(.2,.83,.93,.93,2)
c.GetCanvas().SetRightMargin(0.05)
c.GetCanvas().SetTopMargin(0.05)
c.GetCanvas().SetLeftMargin(0.15)
c.GetCanvas().SetBottomMargin(0.15)
c.cd()

args = parse('--outdir', '--proc')

f =  root.TFile.Open(flatdir + '/%s.root'%args.proc)
t = f.Get('events') 
sel = Selector()

is_top = not(args.proc == 'QCD')
cut = 'gen_pt>200 && fabs(gen_eta)<2.5 && clf_IsMatched'
if is_top:
    cut = tAND(cut, 'gen_size<1.44 && fabs(gen_pdgid)==6')

sel.read_tree(t, cut=cut, 
                branches=['clf_MSD', 'clf_Tau32', 'clf_Tau32SD', 'ptweight'])

nbins = 5
#hbase = root.TH2D('','',nbins,0,400,50,0.2,1.2)
hbase = root.TH2D('','',nbins,0.2,1,50,0,400)
hbase.GetYaxis().SetTitle('m_{SD} [GeV]')
hbase.GetXaxis().SetTitle('#tau_{32}')
hbase.GetXaxis().SetTitleOffset(1); hbase.GetXaxis().SetTitleSize(0.06)
hbase.GetYaxis().SetTitleOffset(1); hbase.GetYaxis().SetTitleSize(0.06)

def draw(key, color):
    fields = ['clf_MSD', key]
    fields.reverse()
    h = sel.draw(fields=fields, hbase=hbase, weight='ptweight')
    h.SetMarkerColor(color)
    h.SetMarkerSize(1)
    h.SetFillColorAlpha(color, 0.1)
    h.SetLineColor(color)
    h.SetMarkerStyle(4)
    return h

h = {
        'Ungroomed' : draw('clf_Tau32', root.kBlue),
        'Groomed' : draw('clf_Tau32SD', root.kRed),
    }


#style = 'VIOLINX(03000130)'
style = 'CANDLEX(10012131)'
#style = 'violin'
h['Ungroomed'].Draw(style)
h['Groomed'].Draw(style+'SAME')
for l,hh in h.iteritems():
    c.GetLegend().AddEntry(hh, l, 'lfp')
c.GetCanvas().Update()

c.Draw(args.outdir, 'msdtau_'+args.proc)

