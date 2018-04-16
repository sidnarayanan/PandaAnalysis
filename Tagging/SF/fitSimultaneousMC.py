#!/usr/bin/env python

import argparse
from sys import argv

NBINS=20

parser = argparse.ArgumentParser(description='fit stuff')
parser.add_argument('--indir',metavar='indir',type=str)
parser.add_argument('--syst',metavar='syst',type=str,default='cent')
args = parser.parse_args()
basedir = args.indir
argv=[]

from math import sqrt
import ROOT as root
from PandaCore.Tools.Load import *
Load('HistogramDrawer')

def imp(w_):
  return getattr(w_,'import')

plot_labels = {
    'smeared' : 'JER smeared',
    'scaleUp' : 'JES Up',
    'scaleDown' : 'JES Down',
    'smearedSJ' : 'SJ JER smeared',
    'scaleSJUp' : 'SJ JES Up',
    'scaleSJDown' : 'SJ JES Down',
    'topPt' : 'No top p_{T} weight',
    'btagUp' : 'b-tag Up',
    'btagDown' : 'b-tag Down',
    'mistagUp' : 'b-mistag Up',
    'mistagDown' : 'b-mistag Down',
    'mergeUp' : 'Merging radius Up',
    'mergeDown' : 'Merging radius Down',
    'cent' : 'Nominal'
    }

plot_postfix = ''
plot_label = plot_labels[args.syst]
if args.syst!='cent':
  if args.syst=='smeared':
    postfix = 'Smeared'
  elif args.syst=='scaleUp':
    postfix = 'ScaleUp'
  elif args.syst=='scaleDown':
    postfix = 'ScaleDown'
  elif args.syst=='smearedSJ':
    postfix = 'Smeared_sj'
  elif args.syst=='scaleSJUp':
    postfix = 'ScaleUp_sj'
  elif args.syst=='scaleSJDown':
    postfix = 'ScaleDown_sj'
  else:
    postfix = ''
  if postfix!='':
    plot_postfix = postfix
  else:
    plot_postfix = args.syst
  args.syst += '_'
else:
  postfix = ''
  args.syst = ''

plot = {}
for iC in [0,1]:
  plot[iC] = root.HistogramDrawer()
  plot[iC].Ratio(1)
  plot[iC].FixRatio(.4)
  plot[iC].DrawMCErrors(False)
  plot[iC].Stack(True)
  plot[iC].DrawEmpty(True)
  plot[iC].SetTDRStyle()
  plot[iC].AddCMSLabel()
  plot[iC].SetLumi(36.6)
  plot[iC].AddLumiLabel(True)
  plot[iC].InitLegend()

hprong = {}; dhprong = {}; pdfprong = {}; norm = {}; smeared = {}; smear = {}; mu = {}; hdata = {}; dh_data={}
mcnorms = {}; mcerrs = {}
mass = root.RooRealVar("m","m_{SD} [GeV]",50,450)

ftemplate = {
      'pass' : root.TFile(basedir+'tag_pass_%shists.root'%args.syst),
      'fail' : root.TFile(basedir+'tag_fail_%shists.root'%args.syst),
    }

inputs = {
    1 : ['QCD','Z+g','Z+q','W+g','W+q','Top [unmatched]'],
    2 : ['Diboson', 'Top [W-matched]'],
    3 : ['Top [t-matched]'],
    }

def build_sum_hist(cat,nprong):
  f = ftemplate[cat]
  hsum = None
  for i in inputs[nprong]:
    if not hsum:
      hsum = f.Get('h_fj1MSD%s_%s'%(postfix,i)).Clone('h_%s_%i'%(cat,nprong))
    else:
      hsum.Add(f.Get('h_fj1MSD%s_%s'%(postfix,i)))
  return hsum

# get histograms
hdata[1] = ftemplate['pass'].Get('h_fj1MSD%s_Data'%postfix)
dh_data[1] = root.RooDataHist('dh_data1','dh_data1',root.RooArgList(mass),hdata[1])
hprong[(3,1)] = build_sum_hist('pass',3) 
hprong[(2,1)] = build_sum_hist('pass',2) 
hprong[(1,1)] = build_sum_hist('pass',1) 

hdata[0] = ftemplate['fail'].Get('h_fj1MSD%s_Data'%postfix)
dh_data[0] = root.RooDataHist('dh_data0','dh_data0',root.RooArgList(mass),hdata[0])
hprong[(3,0)] = build_sum_hist('fail',3) 
hprong[(2,0)] = build_sum_hist('fail',2) 
hprong[(1,0)] = build_sum_hist('fail',1) 

# build pdfs
for iC in [0,1]:
  for iP in xrange(1,4):
    cat = (iP,iC)
    dhprong[cat] = root.RooDataHist('dh%i%i'%cat,'dh%i%i'%cat,root.RooArgList(mass),hprong[cat]) 
    pdfprong[cat] = root.RooHistPdf('pdf%i%i'%cat,'pdf%i%i'%cat,root.RooArgSet(mass),dhprong[cat]) 
    norm_ = hprong[cat].Integral()
    norm[cat] = root.RooRealVar('norm%i%i'%cat,'norm%i%i'%cat,norm_,0.01*norm_,100*norm_)
    mcnorms[cat] = norm_
    err_ = 0
    for iB in xrange(1,hprong[cat].GetNbinsX()+1):
      err_ += pow(hprong[cat].GetBinError(iB),2)
    mcerrs[cat] = sqrt(err_)

# smear pdfs with a single gaussian for pass/fail
jesr_sigma = root.RooRealVar('sigma','sigma',0.1,0.1,10)
jesr_mu = root.RooRealVar('mu','mu',0,-10,10)
#jesr_sigma = root.RooRealVar('jesr_sigma','jesr_sigma',0.1,0.1,0.1)
#jesr_mu = root.RooRealVar('jesr_mu','jesr_mu',0,0,0)
jesr = root.RooGaussian('jesr','jesr',mass,jesr_mu,jesr_sigma)
for iP in [3,2,1]:
  for iC in [0,1]:
    cat = (iP,iC)
    if iP>1 and False:
      smeared[cat] = root.RooFFTConvPdf('conv%i%i'%cat,'conv%i%i'%cat,mass,pdfprong[cat],jesr)
    else:
      smeared[cat] = pdfprong[cat]
'''
for iP in [3,2,1]:
  mu[iP] = root.RooRealVar('mu%i'%iP,'mu%i'%iP,10,-25,25)
  smear[iP] = root.RooGaussian('gauss%i'%iP,'gauss%i'%iP,mass,mu[iP],sigma)
  for iC in [0,1]:
    cat = (iP,iC)
    if iP<3:
      smeared[cat] = pdfprong[cat]
    else:
      smeared[cat] = root.RooFFTConvPdf('conv%i%i'%cat,'conv%i%i'%cat,mass,pdfprong[cat],smear[iP])
'''

model = {}
# nsigtotal = root.RooFormulaVar('nsigtotal','norm30+norm31',root.RooArgList(norm[(3,0)],norm[(3,1)]))
nsigtotal_ = mcnorms[(3,0)]+mcnorms[(3,1)]
nsigtotal = root.RooRealVar('nsigtotal','nsigtotal',nsigtotal_,0.5*nsigtotal_,2*nsigtotal_)
# eff_ = norm[(3,1)].getVal()/(norm[(3,1)].getVal()+norm[(3,0)].getVal())

def calcEffAndErr(p,perr,f,ferr):
  eff_ = p/(p+f)
  err_ = pow( perr * f / pow(p+f,2) , 2 )
  err_ += pow( ferr * p / pow(p+f,2) , 2 )
  err_ = sqrt(err_)
  return eff_,err_

eff_,err_ = calcEffAndErr(mcnorms[(3,1)],mcerrs[(3,1)],mcnorms[(3,0)],mcnorms[(3,1)])
eff_ = mcnorms[(3,1)]/(mcnorms[(3,1)]+mcnorms[(3,0)])
eff = root.RooRealVar('eff','eff',eff_,0.5*eff_,2*eff_)
normsig = {
      0 : root.RooFormulaVar('nsigfail','(1.0-eff)*nsigtotal',root.RooArgList(eff,nsigtotal)),
      1 : root.RooFormulaVar('nsigpass','eff*nsigtotal',root.RooArgList(eff,nsigtotal)),
    }
for iC in [0,1]:
  model[iC] = root.RooAddPdf('model%i'%iC,'model%i'%iC,
#                            root.RooArgList(*[pdfprong[(x,iC)] for x in [1,2,3]]),
                            root.RooArgList(*[smeared[(x,iC)] for x in [1,2,3]]),
                            root.RooArgList(norm[(1,iC)],norm[(2,iC)],normsig[iC]))

# build simultaneous fit
sample = root.RooCategory('sample','')
sample.defineType('pass',1)
sample.defineType('fail',2)

datacomb = root.RooDataHist('datacomb','datacomb',root.RooArgList(mass),
                            root.RooFit.Index(sample),
                            root.RooFit.Import('pass',dh_data[1]),
                            root.RooFit.Import('fail',dh_data[0]))

simult = root.RooSimultaneous('simult','simult',sample)
simult.addPdf(model[1],'pass')
simult.addPdf(model[0],'fail')

# fit!
fitresult = simult.fitTo(datacomb,
                          root.RooFit.Extended(),
                          root.RooFit.Strategy(2),
                          root.RooFit.Minos(root.RooArgSet(eff)),
                          root.RooFit.NumCPU(4),
                          root.RooFit.Save())


# dump the efficiencies
pcat=(3,1); fcat=(3,0)
masslo=110; masshi =210 # to get on binedges for now
effMass_ = hprong[pcat].Integral(hprong[pcat].FindBin(masslo),hprong[pcat].FindBin(masshi)-1)/hprong[pcat].Integral() # prefit efficiency of the mass cut on the pass distribution
mass.setRange('MASSWINDOW',masslo,masshi)
massint = smeared[(3,1)].createIntegral(root.RooArgSet(mass),root.RooArgSet(mass),'MASSWINDOW')
effMass = massint.getVal(); effMassErr = massint.getPropagatedError(fitresult)

# make nice plots

colors = {
  1:8,
  2:6,
  3:root.kOrange
}
labels = {
  1:'Unmatched',
  2:'Matched W',
  3:'Matched top',
}

eff_val = eff.getVal()
nsigtotal_val = nsigtotal.getVal()

for iC in [0,1]:
  for iP in [3,2,1]:
    cat = (iP,iC)
    h = smeared[cat].createHistogram('h%i%i'%cat,mass,root.RooFit.Binning(NBINS))
    h.SetLineWidth(3)
    h.SetLineStyle(1)
    h.SetLineColor(colors[iP])
    if iP!=3:
      h.Scale(norm[cat].getVal()/h.Integral())
    else:
      if iC==0:
        scale = (1-eff_val)*nsigtotal_val
      else:
        scale = eff_val*nsigtotal_val
      h.Scale(scale/h.Integral())

    plot[iC].AddAdditional(h,'hist',labels[iP])
    
  hprefit = hprong[(1,0)].Clone('prefit')
  hprefit.Reset()
  for jP in [3,2,1]:
    hprefit.Add(hprong[(jP,iC)])
  hprefit.SetLineWidth(2)
  hprefit.SetLineStyle(2)
  hprefit.SetLineColor(root.kBlue+2)
  hprefit.SetFillStyle(0)
  plot[iC].AddAdditional(hprefit,'hist','Pre-fit')

  for iP in [1,2,3]:
    cat = (iP,iC)
    hprong[cat].SetLineColor(colors[iP])
    hprong[cat].SetFillStyle(0)
    hprong[cat].SetLineWidth(2)
    hprong[cat].SetLineStyle(2)
    plot[iC].AddAdditional(hprong[cat],'hist')

  hdata[iC].SetLineColor(root.kBlack)
  hdata[iC].SetMarkerStyle(20);
  plot[iC].AddHistogram(hdata[iC],'Data',root.kData)

  hmodel_ = model[iC].createHistogram('hmodel%i'%iC,mass,root.RooFit.Binning(NBINS))
  hmodel = root.TH1D(); hmodel_.Copy(hmodel)
  hmodel.SetLineWidth(3);
  hmodel.SetLineColor(root.kBlue+10)
  hmodel.GetXaxis().SetTitle('fatjet m_{SD} %s [GeV]'%postfix)
  hmodel.GetYaxis().SetTitle('Events/10 GeV')

  signorm_val = (eff_val if iC==1 else 1-eff_val)*nsigtotal_val
  hmodel.Scale((sum([norm[(x,iC)].getVal() for x in [1,2]])+signorm_val)/hmodel.Integral())
  # hmodel.Scale(sum([norm[(x,iC)].getVal() for x in [1,2,3]])/hmodel.Integral())
  hmodel.SetFillStyle(0)
  plot[iC].AddHistogram(hmodel,'Post-fit',root.kExtra5)
  plot[iC].AddAdditional(hmodel,'hist')


  plot[iC].AddPlotLabel('#varepsilon_{tag}^{Data} = %.4g^{+%.2g}_{-%.2g}'%(eff.getVal(),abs(eff.getErrorHi()),abs(eff.getErrorLo())),
                        .6,.47,False,42,.04)
  #plot[iC].AddPlotLabel('#varepsilon_{tag}^{MC} = %.3g^{+%.2g}_{-%.2g}'%(eff_,err_,err_),
  plot[iC].AddPlotLabel('#varepsilon_{tag}^{MC} = %.4g'%(eff_),
                        .6,.37,False,42,.04)
  plot[iC].AddPlotLabel('#varepsilon_{tag+mSD}^{Data} = %.4g^{+%.2g}_{-%.2g}'%(effMass*eff.getVal(),effMass*abs(eff.getErrorHi()),effMass*abs(eff.getErrorLo())),
                        .6,.27,False,42,.04)
  #plot[iC].AddPlotLabel('#varepsilon_{tag+mSD}^{MC} = %.3g^{+%.2g}_{-%.2g}'%(eff_*effMass_,err_*effMass_,err_*effMass_),
  plot[iC].AddPlotLabel('#varepsilon_{tag+mSD}^{MC} = %.4g'%(eff_*effMass_),
                        .6,.17,False,42,.04)

plot[1].AddPlotLabel('Pass category',.18,.77,False,42,.05)
plot[0].AddPlotLabel('Fail category',.18,.77,False,42,.05)
for _,p in plot.iteritems():
  p.AddPlotLabel(plot_label,.18,.7,False,42,.05)
plot_postfix = '_'+plot_postfix
plot[1].Draw(basedir,'/fits/pass%s'%plot_postfix)
plot[0].Draw(basedir,'/fits/fail%s'%plot_postfix)

# save outpuat
w = root.RooWorkspace('w','workspace')
w.imp = imp(w)
w.imp(mass)
for x in [nsigtotal,eff,jesr,sample,datacomb,simult]:
  w.imp(x)
for iC in [0,1]:
  w.imp(normsig[iC])
  w.imp(dh_data[iC])
  w.imp(model[iC])
  for iP in [1,2,3]:
    cat = (iP,iC)
#    w.imp(smear[iP]); w.imp(mu[iP])
    w.imp(dhprong[cat])
    w.imp(pdfprong[cat])
    w.imp(norm[cat])
#    w.imp(smeared[cat])
w.writeToFile(basedir+'/fits/wspace%s.root'%plot_postfix)

s = []
s.append( 'Tagging cut:' )
s.append( '\tPre-fit efficiency was %f'%(eff_) )
s.append( '\tPost-fit efficiency is %f +%.5g -%.5g'%(eff.getVal(),abs(eff.getErrorHi()),abs(eff.getErrorLo())) )
s.append( 'Tagging+mass cut:' )
s.append( '\tPre-fit mass efficiency was %f'%(effMass_) )
s.append( '\tPost-fit mass efficiency is %f +/-%.5g'%(effMass,effMassErr) )
s.append( '\tPre-fit mass+tag efficiency was %f'%(effMass_*eff_) )
s.append( '\tPost-fit mass+tag efficiency is %f +%.5g -%.5g'%(effMass*eff.getVal(),effMass*abs(eff.getErrorHi()),effMass*abs(eff.getErrorLo())) )
s.append( '\tSF = %f +%f -%f'%(effMass*eff.getVal()/(eff_*effMass_), effMass*abs(eff.getErrorHi())/(eff_*effMass_), effMass*abs(eff.getErrorLo())/(eff_*effMass_)) )

s = '\n'.join(s)

with open(basedir+'/fits/summary%s.txt'%plot_postfix,'w') as fout:
  fout.write(s)
print s
