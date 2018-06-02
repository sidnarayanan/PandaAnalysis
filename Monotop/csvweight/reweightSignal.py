#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse

### SET GLOBAL VARIABLES ###
baseDir = getenv('PANDA_FLATDIR')+'/' 
parser = argparse.ArgumentParser(description='reweight CSV shape for fakes')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
args = parser.parse_args()
lumi = 36.56
blind=True
linear=False
sname = argv[0]

argv=[]
import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Utils.load import *
import PandaCore.Tools.Functions
Load('Drawers','HistogramDrawer')

fin = root.TFile(args.outdir+'/ttbar_hists.root')

hdata = fin.Get('h_isojet1CSV_Data')
htoplf = fin.Get('h_isojet1CSV_Top (LF)')
htophf = fin.Get('h_isojet1CSV_Top (HF)')
hqcd = fin.Get('h_isojet1CSV_QCD')
hwj = fin.Get('h_isojet1CSV_W+jets')
hzj = fin.Get('h_isojet1CSV_Z+jets')
hvv = fin.Get('h_isojet1CSV_Diboson')
hfake = htoplf.Clone('h_isojet1CSV_bkg')
hfake.Add(hqcd)
hfake.Add(hwj)
hfake.Add(hzj)
hfake.Add(hvv)
htophf.SetFillStyle(0)

hdata.Add(hfake,-1) # this is what we "estimate" to be real bs
hdata.Scale(1./hdata.Integral())
htophf.Scale(1./htophf.Integral())

hratio = hdata.Clone('hratio')
hratio.Divide(htophf)
hratio.GetYaxis().SetTitle('Data/MC')
hratio.SetMaximum(1.5)
hratio.SetMinimum(0.5)
hratioerr = hratio.Clone()
hratioerr.SetFillStyle(3005)
hratioerr.SetFillColorAlpha(root.kBlack,0.5)
hratioerr.SetLineWidth(0)
hratioerr.SetMarkerStyle(0)

### LOAD PLOTTING UTILITY ###
plot = root.HistogramDrawer()
plot.DrawEmpty(True)
plot.SetLumi(lumi)
plot.SetTDRStyle()
plot.AddCMSLabel()
plot.AddLumiLabel()
plot.SetNormFactor(True)
plot.Logy(True)
plot.SetTDRStyle()
plot.InitLegend()

plotr = root.HistogramDrawer()
plotr.SetLumi(lumi)
plotr.SetRatioStyle()
plotr.AddCMSLabel()
plotr.AddLumiLabel()
plotr.SetGrid()
#plotr.InitLegend(.15,.6,.5,.8)

plot.AddHistogram(htophf,'HF MC',root.kTTbar)
plot.AddHistogram(hdata,'Data-(LF MC)',root.kData)
plot.Draw(args.outdir,'csvweight_tag_shapes')

plotr.AddHistogram(hratio,'Data/MC',root.kExtra1,root.kBlack)
plotr.AddAdditional(hratioerr,'e2')
plotr.Draw(args.outdir,'csvweight_tag_ratio')


fout = root.TFile(args.outdir+'/csvweight_tag.root','RECREATE')
fout.WriteTObject(hdata)
fout.WriteTObject(htophf)
fout.WriteTObject(hratio)
fout.Close()
