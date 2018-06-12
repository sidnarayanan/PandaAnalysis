#!/usr/bin/env python 

import argparse
from sys import argv
from os import getenv
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
parser.add_argument('--process',metavar='process',type=str,default=None)
parser.add_argument('--region',metavar='region',type=str)
args = parser.parse_args()
basedir = getenv('PANDA_FLATDIR')+'/'

argv = []
import ROOT as root
from math import sqrt
from collections import namedtuple 
from array import array
from PandaCore.Utils.load import Load
from PandaCore.Tools.root_interface import draw_hist, read_tree
import PandaAnalysis.VBF.PandaSelection as sel

cut = sel.cuts[args.region]
#'''
for r in ['jot12Mass','fabs(jot12DPhi)','fabs(jot12DEta)']:
    cut = sel.removeCut(cut,r)
#'''
weight = sel.weights[args.region]%35800

Load('PandaCoreTools')
Load('PandaCoreDrawers')

# create some global variables
plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.Logy(True)
plot.InitLegend()
plot.SetLumi(35.8)

plotr = root.HistogramDrawer()
plotr.SetRatioStyle()
plotr.InitLegend(.7,.65,.88,.9)
plotr.SetLumi(35.8)
plotr.SetAutoRange(False)

to_plot = [('1',0,2),('jot1Pt',80,1000),('jot2Pt',40,1000),('jot1Eta',-4.5,4.5),('jot2Eta',-4.5,4.5),
           ('pfmet',200,1000),('jot12Mass',500,5000),('fabs(jot12DPhi)',0,3.14),('fabs(jot12DEta)',0,10)]

repl = { 'nominal':{} }

repl['up'] = {
    'pfmet' : 'pfmetUp',
    'jot1Pt' : 'jot1PtUp',
    'jot2Pt' : 'jot2PtUp',
    'jot1Eta' : 'jot1EtaUp',
    'jot2Eta' : 'jot2EtaUp',
    'jot12Mass' : 'jot12MassUp',
    'jot12DEta' : 'jot12DEtaUp',
    'jot12DPhi' : 'jot12DPhiUp',
}
repl['down'] = {
    'pfmet' : 'pfmetDown',
    'jot1Pt' : 'jot1PtDown',
    'jot2Pt' : 'jot2PtDown',
    'jot1Eta' : 'jot1EtaDown',
    'jot2Eta' : 'jot2EtaDown',
    'jot12Mass' : 'jot12MassDown',
    'jot12DEta' : 'jot12DEtaDown',
    'jot12DPhi' : 'jot12DPhiDown',
}

labels = {
        'WJets' : 'W',
        'ZJets' : 'Z',
        'ZtoNuNu' : 'Z',
        'Diboson' : 'VV'
        }
region_labels = {
        'signal' : 'SR',
        'singlemuon' : '#mu',
        'singleelectron' : 'e',
        'dimuon' : '#mu#mu',
        'dielectron' : 'ee',
        }

f_in = root.TFile.Open(basedir+args.process+'.root')
t_in = f_in.Get('events')

hbase = {}
for d in to_plot:
    nbins = 1 if d[0]=='1' else 20
    hbase[d[0]] = root.TH1D(d[0],'',nbins,d[1],d[2])
    hbase[d[0]].GetYaxis().SetTitle('Events/bin')
    hbase[d[0]].GetXaxis().SetTitle(d[0])


def draw(s):
    weight_ = weight
    cut_ = cut
    to_plot_ = [x[0] for x in to_plot]
    for k,v in repl[s].iteritems():
        cut_ = cut_.replace(k,v)
        for i in xrange(len(to_plot_)):
            to_plot_[i] = to_plot_[i].replace(k,v)
    cut_ = cut_.replace('dphipfmetUp','dphipfmet')
    cut_ = cut_.replace('dphipfmetDown','dphipfmet')
    xarr = read_tree(tree=t_in,branches=to_plot_+[weight_],cut=cut_)
    ret = {}
    for name,shifted_name in zip(to_plot,to_plot_):
        h = hbase[name[0]].Clone(shifted_name)
        draw_hist(h,xarr,[shifted_name],weight_)
        ret[name[0]] = h
    return ret

hists_nominal = draw('nominal')
hists_up = draw('up')
hists_down = draw('down')

print hists_nominal['1'].Integral()
print hists_up['1'].Integral()
print hists_down['1'].Integral()

def plot_yieldo():
    for d in to_plot:
        label_ = labels[args.process]
        label_ += ' (%s)'%(region_labels[args.region])
        plot.cd()
        plot.Reset()
        plot.AddCMSLabel()
        plot.AddLumiLabel(True)
        plot.AddHistogram(hists_nominal[d[0]],'Nominal (stat.)',root.kData)
        plot.AddHistogram(hists_up[d[0]],'JEC Up',root.kExtra2,root.kRed,'hist')
        plot.AddHistogram(hists_down[d[0]],'JEC Down',root.kExtra3,root.kBlue,'hist')
        plot.AddPlotLabel(label_,.18,.77,False,42,.05)
        plotname = 'distribution_%s_%s_%s'%(args.region,args.process,d[0].replace('(','').replace(')',''))
        plot.Draw(args.outdir+'/',plotname)

def plot_ratio():
    for d in to_plot:
        label_ = labels[args.process]
        label_ += ' (%s)'%(region_labels[args.region])
        nominal = hists_nominal[d[0]].Clone()
        up = hists_up[d[0]].Clone()
        down = hists_down[d[0]].Clone()
        up.Divide(nominal)
        down.Divide(nominal)
        for ib in xrange(1,nominal.GetNbinsX()+1):
            val = nominal.GetBinContent(ib)
            nominal.SetBinContent(ib,1)
            if not val: continue
            nominal.SetBinError(ib,nominal.GetBinError(ib)/val)
        for h in [nominal,up,down]:
            h.GetYaxis().SetTitle('Uncertainty')
            h.SetMinimum(0)
            h.SetMaximum(3)
        plotr.cd()
        plotr.Reset()
        plotr.AddCMSLabel()
        plotr.AddLumiLabel(True)
        plotr.AddHistogram(nominal,'Nominal (stat.)',root.kData)
        plotr.AddHistogram(up,'JEC Up',root.kExtra2,root.kRed,'hist')
        plotr.AddHistogram(down,'JEC Down',root.kExtra3,root.kBlue,'hist')
        plotr.AddPlotLabel(label_,.18,.77,False,42,.05)
        plotname = 'distribution_ratio_%s_%s_%s'%(args.region,args.process,d[0].replace('(','').replace(')',''))
        plotr.Draw(args.outdir+'/',plotname)

plot_yieldo()
plot_ratio()
