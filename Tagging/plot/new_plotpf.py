#!/usr/bin/env python

from os import system,getenv
from sys import argv
import argparse
from collections import namedtuple

### SET GLOBAL VARIABLES ###
basedir = getenv('PANDA_FLATDIR')+'/' 
parser = argparse.ArgumentParser(description='plot stuff')
parser.add_argument('--basedir',metavar='basedir',type=str,default=None)
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--cut',metavar='cut',type=str,default='1==1')
parser.add_argument('--region',metavar='region',type=str,default=None)
parser.add_argument('--tt',metavar='tt',type=str,default='')
parser.add_argument('--norm',action='store_true')
args = parser.parse_args()

if args.norm:
    ylabel = 'a.u.'
else:
    ylabel = 'Events/bin'

lumi = 36560.
region = args.region
sname = argv[0]
if args.basedir:
    basedir = args.basedir

argv=[]
import ROOT as root
root.gROOT.SetBatch()
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
from PandaCore.Drawers.plot_utility import *
import PandaAnalysis.Tagging.Selection as sel

### DEFINE REGIONS ###
cut = tAND(sel.cuts[args.region],args.cut)

### LOAD PLOTTING UTILITY ###
plot = PlotUtility()
plot.Stack(True)
plot.Ratio(True)
plot.FixRatio(0.4)
plot.SetTDRStyle()
plot.InitLegend()
plot.DrawMCErrors(True)
plot.AddCMSLabel()
plot.cut = cut
plot.SetLumi(lumi/1000)
plot.SetNormFactor(args.norm)
plot.AddLumiLabel(True)
plot.do_overflow = True
plot.do_underflow = True
plot.AddPlotLabel("m_{SD} > 50 GeV",.18,.77,False,42,.04)

weight = sel.weights[region]%lumi
plot.mc_weight = weight

### DEFINE PROCESSES ###
wjetsq         = Process('W+q',root.kWjets); wjetsq.additional_cut = 'abs(fj1HighestPtGen)!=21'
wjetsg         = Process('W+g',root.kExtra2); wjetsg.additional_cut = 'abs(fj1HighestPtGen)==21'

zjetsq         = Process('Z+q',root.kZjets); zjetsq.additional_cut = 'abs(fj1HighestPtGen)!=21'
zjetsg         = Process('Z+g',root.kExtra5); zjetsg.additional_cut = 'abs(fj1HighestPtGen)==21'

gjetsq         = Process('#gamma+q',root.kGjets); gjetsq.additional_cut = 'abs(fj1HighestPtGen)!=21'
gjetsg         = Process('#gamma+g',root.kExtra3); gjetsg.additional_cut = 'abs(fj1HighestPtGen)==21'

diboson        = Process('Diboson',root.kDiboson)
qcd             = Process("QCD",root.kQCD)

top            = Process('Top [t-matched]',root.kTTbar)
top.additional_cut = 'fj1IsMatched==1 && fj1GenSize<1.44'

wtop           = Process('Top [W-matched]',root.kExtra1)
wtop.additional_cut = 'fj1IsMatched==1 && fj1GenSize>1.44 && fj1IsWMatched==1 && fj1GenWSize<1.44'

untop           = Process('Top [unmatched]',root.kExtra4)
untop.additional_cut = '!((fj1IsMatched==1&&fj1GenSize>1.44 && fj1IsWMatched && fj1GenWSize<1.44)||(fj1IsMatched==1&&fj1GenSize<1.44))'

data            = Process("Data",root.kData)
if region=='photon':
    processes = [qcd,gjetsg,gjetsq,data]
elif region=='mistag':
    processes = [qcd,diboson,zjetsg,zjetsq,untop,wtop,top,wjetsg,wjetsq,data]
elif region=='tag':
    processes = [qcd,diboson,zjetsg,zjetsq,wjetsg,wjetsq,untop,wtop,top,data]
else:
    processes = [qcd,diboson,wjetsg,wjetsq,untop,wtop,top,zjetsg,zjetsq,data]

### ASSIGN FILES TO PROCESSES ###
for p in [wjetsq,wjetsg]:
    p.add_file(basedir+'WJets.root')

for p in [zjetsq,zjetsg]:
    p.add_file(basedir+'ZJets.root')

for p in [gjetsq,gjetsg]:
    p.add_file(basedir+'GJets.root')

for p in [top,wtop,untop]:
    p.add_file(basedir+'TTbar.root')
    p.add_file(basedir+'SingleTop.root')

diboson.add_file(basedir+'Diboson.root')

if 'pho' in region:
    qcd.add_file(basedir+'SinglePhoton.root')
    qcd.additional_cut = sel.triggers['pho']
    qcd.use_common_weight = False
    qcd.additional_weight = 'sf_phoPurity'
    data.add_file(basedir+'SinglePhoton.root')
    data.additional_cut = sel.triggers['pho']
else:
    qcd.add_file(basedir+'QCD.root')
    data.add_file(basedir+'MET.root')
    data.additional_cut = sel.triggers['met']

for p in processes:
    plot.add_process(p)

plot.add_distribution(FDistribution('fj1MSD',50,550,20,'fatjet m_{SD} [GeV]',ylabel))
plot.add_distribution(FDistribution('fj1Pt',200,1000,20,'fatjet p_{T} [GeV]',ylabel))
plot.add_distribution(FDistribution('top_ecf_bdt',-1,1,20,'Top BDT',ylabel))
#plot.add_distribution(FDistribution('top_ecfAll_bdt',-1,1,20,'Top 50ECF BDT',ylabel))
plot.add_distribution(FDistribution('fj1MaxCSV',0,1,20,'fatjet max CSV',ylabel))
plot.add_distribution(FDistribution("1",0,2,1,"dummy","dummy"))

class Ratio:
    def __init__(self,num,den,x,lo,hi,idx):
        self.num=num
        self.den=den
        self.x=x
        self.lo=lo
        self.hi=hi
        self.idx=idx
    def convert(self,t):
        return (t[0],t[1],int(10*t[2]))
    def formula(self):
        snum = 'fj1ECFN_%i_%i_%.2i'%self.convert(self.num)
        sden = 'fj1ECFN_%i_%i_%.2i'%self.convert(self.den)
        return '%s/pow(%s,%.3f)'%(snum,sden,self.x)
    def label(self):
        snum = 'e(%i,%i,%.1f)'%self.num
        sden = 'e(%i,%i,%.1f)'%self.den
        if self.x==1:
            return '%s/%s'%(snum,sden)
        else:
            return '%s/%s^{%.2f}'%(snum,sden,self.x)
    def filename(self):
        return 'input%i'%(self.idx)

ratios = [
            # these are misc
            Ratio((1,3,2),(2,3,1),1,0,1,50),
            Ratio((2,4,.5),(1,3,.5),2,1,3,274),
            Ratio((1,2,4),(1,2,2),2,3,12,11),
            Ratio((1,3,4),(1,2,1),4,0,5,58),
            Ratio((1,3,1),(1,3,.5),2,7,25,31),
            Ratio((1,4,2),(1,4,4),.5,0,.2,250),

            # these are BDT inputs
            Ratio((1,2,2),(1,2,1),2,2,10,0),
            Ratio((1,3,4),(2,3,2),1,0,1,1),
            Ratio((3,3,1),(1,3,4),0.75,0.5,5,2),
            Ratio((3,3,1),(2,3,2),.75,0.4,1.4,3),
            Ratio((3,3,2),(3,3,4),.5,0,.25,4),
            Ratio((1,4,2),(1,3,1),2,0,2,5),
            Ratio((1,4,4),(1,3,2),2,0,2.5,6),
            Ratio((2,4,.5),(1,3,.5),2,1.3,2.5,7),
            Ratio((2,4,1),(1,3,1),2,1,4,8),
            Ratio((2,4,1),(2,3,.5),2,0,1.52,9),
            Ratio((2,4,2),(1,3,2),2,0,5,10),
        ]

#for r in ratios:
#    plot.add_distribution(FDistribution(r.formula(),r.lo,r.hi,20,r.label(),ylabel,filename=r.filename()))

### DRAW AND CATALOGUE ###
plot.draw_all(args.outdir+'/'+region+'_')
