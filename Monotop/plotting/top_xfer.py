#!/usr/bin/env python

from os import system,getenv
from sys import argv
from array import array
import argparse

### SET GLOBAL VARIABLES ###
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cat',metavar='cat',type=str)
args = parser.parse_args()
lumi = 35800.
blind=True
linear=False
sname = argv[0]

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Utils.load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
import PandaAnalysis.Monotop.CombinedSelection as sel
from PandaCore.Tools.root_interface import draw_hist,read_tree
from PandaCore.Drawers.plot_utility import divide_bin_width
Load('HistogramDrawer')

add_cut = 'top_ecf_bdt<0.45' if args.cat=='loose' else 'top_ecf_bdt>0.45'

recoil_bins = array('f',[250,280,310,350,400,450,600,1000])
hbase = root.TH1D('hbase','',len(recoil_bins)-1,recoil_bins)
hbase.GetXaxis().SetTitle('U [GeV]')
hbase.GetYaxis().SetTitle('(Alternate)/(Nominal)')
hbase.SetMinimum(0.95)
hbase.SetMaximum(1.1)

plotr = root.HistogramDrawer()
plotr.SetRatioStyle()
plotr.InitLegend()
plotr.AddCMSLabel()
plotr.SetLumi(lumi/1000)
plotr.AddLumiLabel(True)
plotr.SetAutoRange(False)

fin = root.TFile(getenv('PANDA_FLATDIR')+'/TTbar.root')
tree = fin.Get('events')

def divide(hnum,hden):
    hratio = hbase.Clone()
    for ib in xrange(1,hratio.GetNbinsX()+1):
        v_den = hden.GetBinContent(ib)
        e_num = hnum.GetBinError(ib)
        v_num = hnum.GetBinContent(ib)
        hratio.SetBinContent(ib,v_num/v_den)
        hratio.SetBinError(ib,e_num/v_den)
    return hratio

def build_hists(cut,weight,observable):
    h = {}
    weight = weight%1
    cut = tAND(cut,add_cut)
    weights = {
            '13TeV' : weight,
            '8TeV' : weight.replace('sf_tt','sf_tt8TeV'),
            'none' : weight.replace('sf_tt','1')
    }
    xarr = read_tree(tree,[observable]+weights.values(),cut)
    for label,weight in weights.iteritems():
        h[label] = hbase.Clone(label+'_'+observable)
        draw_hist(h[label],xarr,[observable],weight)
        divide_bin_width(h[label])

    return h



h_sig = build_hists(sel.cuts['signal'],sel.weights['signal'],'pfmet')
h_con = build_hists(sel.cuts['singlemuontop'],sel.weights['singlemuontop'],'pfUWmag')

for k in h_sig:
    h_sig[k].Divide(h_con[k])

h_sig['8TeV'] = divide(h_sig['8TeV'],h_sig['13TeV'])
h_sig['none'] = divide(h_sig['none'],h_sig['13TeV'])
h_sig['13TeV'] = divide(h_sig['13TeV'],h_sig['13TeV'])

plotr.AddHistogram(h_sig['13TeV'],'Top p_{T} weight [nominal]',root.kData,root.kRed)
#plotr.AddHistogram(h_sig['8TeV'],'8 TeV',root.kExtra1)
plotr.AddHistogram(h_sig['none'],'No weight',root.kExtra2)

plotr.Draw(args.outdir,'/top_xfer_pt_'+args.cat)
