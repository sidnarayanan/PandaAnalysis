#!/usr/bin/env python

import argparse
from sys import argv
import json 

parser = argparse.ArgumentParser(description='summary')
parser.add_argument('--indir',metavar='indir',type=str)
parser.add_argument('--wp',type=str,default='tight')
args = parser.parse_args()
basedir = args.indir
argv=[]

from math import sqrt
import ROOT as root
from PandaCore.Utils.load import *
from PandaAnalysis.Tagging.TnPSel import pt_cut, pt_bins
from array import array
Load('HistogramDrawer')

NPT = len(pt_bins) - 1
h_sf = root.TH1D('h','h',NPT,0,1)
h_sf.GetXaxis().SetTitle('p_{T} [GeV]')
h_sf.GetYaxis().SetTitle('#epsilon_{Data}/#epsilon_{MC}')
h_sf.SetLineColor(root.kWhite)
sfs = {x:[] for x in ['cent', 'stathi', 'systhi', 'statlo', 'systlo']}
dump = ['%10s, %6s, %6s, %6s, %6s, %6s'%('ptbin', 'sf', 'stathi', 'statlo', 'systhi', 'systlo')]
zeros = []
xs = []

for i in xrange(NPT):
  h_sf.GetXaxis().SetBinLabel(i+1,'%i-%i'%tuple(pt_bins[i:i+2]))
  data = json.load(open(args.indir+'/summary__%s_%i.json'%(args.wp,i)))
  h_sf.SetBinContent(i+1,data['sf'])
  sfs['cent'].append(data['sf'])
  sfs['stathi'].append(data['hi'])
  sfs['statlo'].append(data['lo'])
  sfs['systhi'].append(data['hi'])
  sfs['systlo'].append(data['lo'])
  zeros.append(0.01)
  xs.append((i+0.5)/NPT)
  dump.append('%10s, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f'%('%i-%i'%tuple(pt_bins[i:i+2]), data['sf'],
                                                         data['hi'], data['lo'], data['hi'],
                                                         data['lo']))
zeros = array('f',zeros)
xs = array('f',xs)
  

g_stat = root.TGraphErrors(NPT,xs,
                           array('f', sfs['cent']), zeros,
                           array('f', sfs['stathi']))
g_syst = root.TGraphErrors(NPT,xs,
                           array('f', sfs['cent']), zeros,
                           array('f', sfs['stathi']))
g_syst.SetLineColor(root.kRed)
g_syst.SetFillColor(root.kRed)

with open(args.indir+'/summary_%s.txt'%args.wp,'w') as fdump:
  fdump.write('\n'.join(dump))

plot = root.HistogramDrawer()
plot.SetTDRStyle()
plot.InitLegend()
root.gStyle.SetOptStat(0)
plot.SetAbsMin(0.8)

plot.AddHistogram(h_sf,'',10,10,'hist')
plot.AddAdditional(g_syst,'2','Sys.+stat.')
plot.AddAdditional(g_stat,'e0','Stat.')
plot.Draw(args.indir,'summary_'+args.wp)
