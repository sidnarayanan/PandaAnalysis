#!/usr/bin/env python
from sys import argv,exit
import argparse
from os import getenv

parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default='.')
parser.add_argument('--cut',metavar='cut',type=str,default=None)
parser.add_argument('--pt',metavar='pt',type=str,default='inc')
args = parser.parse_args()

figsdir = args.outdir
sname = argv[0]

argv = []
import ROOT as root
from PandaCore.Tools.Load import *
cut = '' 
plotlabel = '110 < m_{SD} < 210 GeV'
xcut = 'mSD>110 && mSD<210 && pt<1000'

Load('PandaCoreDrawers')
roc = root.ROCTool()

fin = root.TFile(figsdir+'/'+'hists.root')

roc.Logy()
roc.SetPlotRange(0.005,1)
roc.InitCanvas(.7,.2,.85,.45, False)
roc.SetFile(fin)
roc.c.AddPlotLabel(plotlabel,.2,.82,False,42,.04)

variables = [
#  ('top_ecfv14_bdt','ECF+#tau_{32}^{SD}+f_{rec} BDT v2',1,3),
#  ('top_ecfv12_bdt','ECF BDT v2',2,3),
#  ('top_ecfv8_bdt','ECF+#tau_{32}^{SD}+f_{rec} BDT',1,1),
  ('top_ecfv8_bdt','BDT',1,1),
#  ('top_ecfv6_bdt','ECF BDT',2,1),
#  ('top_taufrec_bdt','#tau_{32}^{SD}+f_{rec} BDT',4,1),
#  ('htt_frec','f_{rec}',1,2),
  ('tau32SD','#tau_{3}/#tau_{2}',3,2),
#  ('tau32','#tau_{32}',4,2),
#  ('input0','e(1,2,2)/e(1,2,1)^{2}',3,3),
#  ('input1','e(1,3,4)/e(2,3,2)',4,3),
#  ('input2','e(3,3,1)/e(1,3,4)^{3/4}',5,3),
#  ('input3','e(3,3,1)/e(2,3,2)^{3/4}',6,3),
#  ('input4','e(3,3,2)/e(3,3,4)^{1/2}',7,3),
#  ('input5','e(1,4,2)/e(1,3,1)^{2}',8,3),
#  ('input6','e(1,4,4)/e(1,3,2)^{2}',9,3),
#  ('input7','e(2,4,0.5)/e(1,3,0.5)^{2}',10,3),
#  ('input8','e(2,4,1)/e(1,3,1)^{2}',11,3),
#  ('input9','e(2,4,1)/e(2,3,0.5)^{2}',12,3),
#  ('input10','e(2,4,2)/e(1,3,2)^{2}',13,3),
  # ('N3_05','N3, #beta=0.5',11,2),
  # ('N3_10','N3, #beta=1.0',6,2),
  # ('N3_20','N3, #beta=2.0',10,2),
            ] 
for iV in xrange(len(variables)):
  PInfo(sname,variables[iV][1])
  if len(variables[iV])==2:
    v,vlabel = variables[iV]
    vcolor = iV
    vstyle = 1 if 'bdt' in v else 2
  else:
    v,vlabel,vcolor,vstyle = variables[iV]
  roc.CalcROC('h_%s_Top quark jets'%v,'h_%s_qOverg jets'%v,vlabel,vcolor,vstyle,1)

base = getenv('SCRAMJETFLAT')+'/'
fsig = root.TFile(base+'ZpTT.root'); tsig = fsig.Get('puppiCA15')
fbkg = root.TFile(base+'QCD_evt8.root'); tbkg = fbkg.Get('puppiCA15')
#roc.CalcWP(tsig,tbkg,1,tAND(xcut,'matched==1&&gensize<1.2'),'ptweight*normalizedWeight',xcut,'ptweight_analytic*normalizedWeight','tau32SD<0.6&&fabs(htt_frec)<0.2','#tau_{32}^{SD}<0.6 && f_{rec}<0.2')

roc.DrawAll(figsdir,cut+'roc')

