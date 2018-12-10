#!/usr/bin/env python

from PandaCore.Tools.script import * 
from re import sub 
from collections import namedtuple
import numpy as np
from sys import exit 

Load('CanvasDrawer')

args = parse('--outdir')

f = root.TFile.Open(args.outdir + '/singlemuonw_hists.root')

keys = [k.GetName() for k in f.GetListOfKeys()]
inputs = sorted([int(k.replace('h_ratio_','').replace('_Data','')) for k in keys
                 if '_Data' in k and '_ratio_' in k])

test = "e(1,2,0.5)/e(1,2,1.0)^{0.50}"
ECFR = namedtuple('ECFR', ['a','N','alpha','b','M','beta'])

def get_params(h):
    if type(h) == str:
        title = h
    else:
        title = h.GetXaxis().GetTitle()
    title = sub(r'\)\)\^.*','',title)
    title = sub(r'e\(','',title)
    title = sub(r'\)/\(',',',title)
    title = sub(r'\)','',title)
    return ECFR(*[float(x.strip()) for x in title.split(',')])

def compute_ks(i):
    h_data = f.Get('h_ratio_%i_Data'%i)
    h_mc = None 
    for k in keys:
        if '_Data' in k:
            continue
        if '_ratio_%i_'%i in k:
            if h_mc is None:
                h_mc = f.Get(k)
            else:
                h_mc.Add(f.Get(k))

    n_nonzero = sum([1 for x in xrange(1, h_data.GetNbinsX()+1)
                       if h_data.GetBinContent(x) > 0])

    if h_data.Integral() == 0 or n_nonzero < 3:
        ks = 0
    else:
        ks = h_data.KolmogorovTest(h_mc)
    return ks, h_data

kss = {}
for i in inputs:
    ks,hdata = compute_ks(i)
    kss[get_params(hdata)] = ks

params = {x:set([]) for x in ECFR._fields}
for e in kss:
    for p in params:
        params[p].add(getattr(e, p))

for p in params:
    params[p] = sorted(params[p])


c = root.CanvasDrawer(800, 600)
c.SetTDRStyle()
c.GetCanvas().SetRightMargin(0.16)
c.GetCanvas().SetLeftMargin(0.1)
c.cd()
for a in params['a']:
    for b in params['b']:
        c.Reset()
        h2 = root.TH2D('','',3, 1.5, 4.5, 3, 1.5, 4.5)
        ks = []
        for N in params['N']:
            for M in params['M']:
                for alpha in params['alpha']:
                    for beta in params['beta']:
                        e = ECFR(a, N, alpha, b, M, beta)
                        if e in kss:
                            found = True
                            ks.append(kss[e])
                if ks:
                    h2.SetBinContent(h2.FindBin(N, M), -np.log10(sum(ks)/len(ks)))
        h2.SetMinimum(0)
#        h2.SetMinimum(-15)
        h2.Draw('colz text')
        c.Draw(args.outdir, 'scanMN_%i_%i'%(int(a), int(b)))


c.Reset()
h2 = root.TH2D('','',3, 1.5, 4.5, 3, 1.5, 4.5)
for N in params['N']:
    for M in params['M']:
        if M > N:
            continue 
        ks = []
        for a in params['a']:
            for b in params['b']:
                for alpha in params['alpha']:
                    for beta in params['beta']:
                        e = ECFR(a, N, alpha, b, M, beta)
                        if e in kss:
                            if kss[e] > 0:
                                ks.append(kss[e])
                            if kss[e] == 1:
                                pass
        if ks:
            h2.SetBinContent(h2.FindBin(N, M), -np.log10(sum(ks)/len(ks)))
#        h2.SetMaximum(0)
#        h2.SetMinimum(-15)
h2.Draw('colz text')
c.Draw(args.outdir, 'scanMN')

N = 4.; M = 3.
c.Reset()
h2 = root.TH2D('','',3, .5, 3.5, 3, .5, 3.5)
for a in params['a']:
    for b in params['b']:
        ks = []
        for alpha in params['alpha']:
            for beta in params['beta']:
                e = ECFR(a, N, alpha, b, M, beta)
                if e in kss:
                    if kss[e] > 0:
                        ks.append(kss[e])
                else:
                    pass
        if ks:
            h2.SetBinContent(h2.FindBin(a, b), -np.log10(sum(ks)/len(ks)))
#h2.SetMaximum(0)
#h2.SetMinimum(-15)
h2.Draw('colz text')
c.Draw(args.outdir, 'scanab')


order =['2/2', '3/3', '4/4', '4/3', '3/2', '4/2']
# x
h2 = root.TH2D('','',6,-0.5,5.5,20,0.5,8)
g = root.TGraph2D()
g1 = root.TGraph()
for i,o in enumerate(order):
    h2.GetXaxis().SetBinLabel(i+1, o)
pts = {}
for e,ks in kss.iteritems():
    if ks <= 0 or ks == 1:
        continue
    o = '%i/%i'%(int(e.N), int(e.M))
    x = e.a * e.alpha / (e.b * e.beta)
    t = (o,x)
    if t not in pts:
        pts[t] = []
    pts[t].append(ks)

i = 0
for p,v in pts.iteritems():
    vv = -np.log10(sum(v)/len(v))
    g.SetPoint(i, order.index(p[0]), p[1], vv)
    g1.SetPoint(i, order.index(p[0]), p[1])
    i += 1

for ix in xrange(1, h2.GetNbinsX()+1):
    for iy in xrange(1, h2.GetNbinsY()+1):
        x = h2.GetXaxis().GetBinCenter(ix)
        y = h2.GetYaxis().GetBinCenter(iy)
        h2.SetBinContent(ix, iy, g.Interpolate(x, y))
h2.Draw('colz')
h2.GetXaxis().SetTitle('N/M')
h2.GetYaxis().SetTitle('a#alpha/b#beta')
h2.GetZaxis().SetTitle('-log_{10}(KS)')
h2.SetMinimum(0)
g1.SetMarkerStyle(20)
#g1.Draw('p same')
c.Draw(args.outdir, 'scanmnab')
