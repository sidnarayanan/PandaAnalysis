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

#determine the SF from low dphi data
def sf(hdata, hqcd, hbkg):
    h = hdata.Clone()
    h.Add(hbkg, -1)
    h.GetYaxis().SetTitle('(Data-Bkg)/QCD')
    h.Divide(hqcd)
    return h 

def sfs():
    r = {}
    hdata = get(flow, 'Data')
    hqcd = get(flow, 'QCD')

    hbkg = None 
    for v in ['W', 'Z']:
        for d in ['EWK','QCD']:
            h_ = get(flow, '%sPlusjets [%s]'%(v,d))
            if not hbkg:
                hbkg = h_.Clone('hbkg')
            else:
                hbkg.Add(h_)

    r['nominal'] = sf(hdata, hqcd, hbkg)

    # shift theory up/down 20%
    hbkg.Scale(0.8)
    r['wz_down'] = sf(hdata, hqcd, hbkg)
    hbkg.Scale(1.5)
    r['wz_up'] = sf(hdata, hqcd, hbkg)

    return r 

def shapes(r):
    s = {} 
    htarget = get(fhigh, 'QCD')
    for k,v in r.iteritems():
        h = htarget.Clone() 
        h.GetYaxis().SetTitle('Events/GeV')
        h.Multiply(v)
        s[k] = h 

    return s 

def interp(r):
    i = {} 
    bins = array('f', [200, 400, 600, 900, 1200, 1500, 2000, 2750, 3500, 5000])
    bin_centers = [(x+y)*0.5 for x,y in zip(bins[1:],bins[:-1])]
    hbase = root.TH1D('hbase','hbase',len(bins)-1, bins)
    for k,v in r.iteritems():
        h = hbase.Clone()
        v_err = v.Clone()
        for ib in xrange(1, v.GetNbinsX()):
            central = v.GetBinContent(ib)
            error = v.GetBinError(ib)
            v_err.SetBinContent(ib, error / central ) # relative error 
        g = root.TGraphErrors(v) # how we interpolate
        g_err = root.TGraphErrors(v_err) 
        for ib, c in enumerate(bin_centers):
            h.SetBinContent(ib + 1, g.Eval(c))
            h.SetBinError(ib + 1, h.GetBinContent(ib + 1) * g_err.Eval(c))
        h.GetYaxis().SetTitle('Events/GeV')
        h.GetXaxis().SetTitle('m_{jj} [GeV]')
        i[k] = h 
    return i 


plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.Logy(False)
plot.InitLegend()
plot.SetLumi(36)

correction = sfs()
plot.AddCMSLabel()
plot.AddLumiLabel(True)
plot.AddHistogram(correction['nominal'], 'Nominal', root.kData, -1, 'hist e')
plot.AddHistogram(correction['wz_up'], 'W/Z theory unc', root.kZjets, -1, 'hist')
plot.AddHistogram(correction['wz_down'], '', root.kZjets, -1, 'hist')
plot.Draw(args.outdir, 'correction')

s = shapes(correction)
plot.Logy(True)
plot.Reset(True)
plot.AddCMSLabel()
plot.AddLumiLabel(True)
plot.AddHistogram(s['nominal'], 'Nominal', root.kData, -1, 'hist e')
plot.AddHistogram(s['wz_up'], 'W/Z theory unc', root.kZjets, -1, 'hist')
plot.AddHistogram(s['wz_down'], '', root.kZjets, -1, 'hist')
plot.Draw(args.outdir, 'template_coarse')


s = interp(s)
plot.Logy(True)
plot.Reset(True)
plot.AddCMSLabel()
plot.AddLumiLabel(True)
plot.AddHistogram(s['nominal'], 'Nominal', root.kData, -1, 'hist e')
plot.AddHistogram(s['wz_up'], 'W/Z theory unc', root.kZjets, -1, 'hist')
plot.AddHistogram(s['wz_down'], '', root.kZjets, -1, 'hist')
plot.Draw(args.outdir, 'template')

