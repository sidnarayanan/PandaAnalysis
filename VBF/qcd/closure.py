#!/usr/bin/env python

from os import system,mkdir,getenv
from sys import argv
from array import array
import argparse
from numpy import arange, sqrt

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str)
args = parser.parse_args()

figsdir = args.outdir
argv=[]

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import Load 
Load('HistogramDrawer')

flow = root.TFile.Open(args.outdir + '/lowdphi_hists.root')
fhigh = root.TFile.Open(args.outdir + '/highdphi_hists.root')

def get(f, p):
    return f.Get('h_jot12Mass_%s'%p)

setattr(flow, 'get', lambda x : get(flow, x))
setattr(fhigh, 'get', lambda x : get(fhigh, x))

#determine the SF from low dphi data
def sf(f):
    f.ls()
    hdata = f.Get('h_jot12Mass_Data')
    print hdata
    print hdata.GetBinContent(1)
    hqcd = f.get('QCD')
    hbkg = None 
    for v in ['W', 'Z']:
        for d in ['EWK','QCD']:
            h_ = f.get('%sPlusjets [%s]'%(v,d))
            if not hbkg:
                hbkg = h_.Clone('hbkg')
            else:
                hbkg.Add(h_)
    print f
    print hdata.GetBinContent(1), hbkg.GetBinContent(1), hqcd.GetBinContent(1)
    h = hdata.Clone()
    h.Add(hbkg, -1)
    h.GetYaxis().SetTitle('(Data-Bkg)/QCD')
    h.Divide(hqcd)
    return h 

    return r 



plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.Logy(False)
plot.InitLegend()
plot.SetLumi(36)

sf_high = sf(fhigh)
sf_low = sf(flow)

print sf_low.GetBinContent(1)
print sf_high.GetBinContent(1)

ratio = sf_low.Clone(); ratio.Divide(sf_high)
ratio.GetXaxis().SetTitle('m_{jj} [GeV]')
ratio.GetYaxis().SetTitle('(SF low #Delta#phi)/(SF high #Delta#phi)')

plot.AddCMSLabel()
plot.AddLumiLabel(True)
plot.AddHistogram(ratio, '', root.kData, -1, 'hist e')
plot.Draw(args.outdir, 'closure')

