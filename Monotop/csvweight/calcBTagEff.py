#!/usr/bin/env python

from sys import argv,exit
import argparse
from array import array

parser = argparse.ArgumentParser(description='skip stuff')
parser.add_argument('--infile',type=str)
parser.add_argument('--outdir',type=str)
parser.add_argument('--wp',type=str)
parser.add_argument('--flav',type=int)
args = parser.parse_args()
sname = argv[0]
argv = []

from PandaCore.Utils.load import *
from PandaCore.Tools.Misc import *
import ROOT as root

Load('Drawers','CanvasDrawer')

plot = root.CanvasDrawer()
plot.SetTDRStyle()
root.gStyle.SetPadRightMargin(0.15)
root.gStyle.SetNumberContours(999);
root.gStyle.SetPalette(root.kBird)
root.gStyle.SetPaintTextFormat(".2g")
c = root.TCanvas()

f = root.TFile(args.infile)
jets = f.Get('jets')

wps = {
		'L' : 0.5426,
		'M' : 0.8484
		}

aeta = array('f',[0,0.5,1.5,2.5]); neta = len(aeta)-1
apt = array('f',[20,50,80,120,200,300,400,500,700,1000]); npt = len(apt)-1

hbase = root.TH2F('hbase','hbase',neta,aeta,npt,apt)
hnum = hbase.Clone('hnum')

hbase.Sumw2()
hnum.Sumw2()

jets.Draw('pt:fabs(eta)>>hbase','flavor==%i'%(args.flav))
jets.Draw('pt:fabs(eta)>>hnum','csv>%f && flavor==%i'%(wps[args.wp],args.flav))

hnum.Divide(hbase)

hnum.GetYaxis().SetTitle('p_{T} [GeV]')
hnum.GetXaxis().SetTitle('|#eta|')
hnum.GetZaxis().SetTitle('#epsilon(CSV>%.3f, flavor=%i)'%(wps[args.wp],args.flav))

c.cd()
hnum.Draw('colz text')
plot.AddCMSLabel(.16,.94)
plot.SetLumi(36.6); plot.AddLumiLabel(True)
plot.SetCanvas(c)

plot.Draw(args.outdir,'/btag_eff_%sWP_flav%i'%(args.wp,args.flav))

etastr = '{' + ','.join(map(str,aeta)) + '}'
ptstr = '{' + ','.join(map(str,apt)) + '}'

print 'eta bins = ',etastr
print 'pt bins = ',ptstr

effstrs = []
for iE in xrange(1,neta+1):
	effs = []
	for iP in xrange(1,npt+1):
		effs.append('%.3f'%hnum.GetBinContent(iE,iP))
	effstrs.append( '{' + ','.join(effs) + '}' )

effstr = '{' + ',\n'.join(effstrs) + '}'
print 'effs = '
print effstr
