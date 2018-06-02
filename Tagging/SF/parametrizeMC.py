#!/usr/bin/env python

from sys import argv,exit
argv=[]

import ROOT as root
from PandaCore.Utils.load import *
import PandaCore.Statistics.RooFitUtils as RooFitUtils 
from PandaCore.Statistics.SimpleStats import SimpleVar

basedir = '/home/snarayan/home000/store/panda/v8/sf/'
nbins = 40

mass_var = SimpleVar('mSD',50,350,title='m_{SD} [GeV]')
mass = mass_var.rvar
mass.setBins(10000,'cache') # for FFT
mass.setMin('cache',50)
mass.setMax('cache',350)

tag = SimpleVar('top_ecfv8_bdt',-1.2,1)
weight = SimpleVar('weight',0,4)

cat = RooFitUtils.CategoryManager('cat')

def fitMC(nprongs,isPass,model='GAUSSEXPERF',smear=None,paramset=None):
  fmc = root.TFile(basedir+'/tag__prong%i.root'%nprongs); tmc = fmc.Get('prong%i'%nprongs)
  if isPass:
    cut = 'top_ecfv8_bdt>0.66'
  else:
    cut = 'top_ecfv8_bdt<0.66'
  ds = RooFitUtils.treeToDS(tmc,[mass,tag.rvar,weight.rvar],cut=cut,weight='weight')

  hmass = root.TH1D('hmass','hmass',nbins,50,350)
  tmc.Draw('mSD>>hmass','weight*(%s)'%cut,'')
  integral = hmass.Integral() # get the integral so we can do an extended likelihood by default

  smear_params = {
      'mu':(0,-20,20),
      'sigma':(10,.1,100),
      'n_cb' : (1,0.01,5),
      'alpha':(-1,-100,0),
      }

  dparams = {}
  dparams['fail'] = {
      'a':(140,-1000,1000),
      'b':(58.1,0.01,1000),
      'c':(-0.0155,-10,0),
      'alpha':(-1.4,-10,10),
      'mu':(165,145,185),
      'sigma':(20,10,100),
      'n_cb' : (1.7,0,5),
      'norm_gauss':(integral,1,integral*2),
      'norm_bw':(integral,1,integral*2),
      'norm_cb':(integral,1,integral*2),
      'norm_experf':(integral*0.3,1,integral*2),
      'w_mu':(80,78,82),
      'w_sigma':(8,5,20),
      }

  dparams['pass'] = {
      'a':(200,-1000,1000),
      'b':(50,0.01,1000),
      'c':(-0.0155,-10,0),
      'alpha':(-1.4,-10,10),
      'mu':(165,145,185),
      'sigma':(20,10,100),
      'n_cb' : (1.7,0,5),
      'norm_gauss':(integral,1,integral*2),
      'norm_bw':(integral,1,integral*2),
      'norm_cb':(integral,1,integral*2),
      'norm_experf':(integral*0.2,1,integral*2),
      'w_mu':(80,78,82),
      'w_sigma':(8,5,20),
      }

  dparams['fail2'] = {
      'a':(200,-1000,1000),
      'b':(5.6,0.01,1000),
      'c':(-0.0155,-10,0),
      'alpha':(-1.4,-10,10),
      'mu':(80,70,90),
      'sigma':(20,10,100),
      'n_cb' : (1.7,0,5),
      'norm_gauss':(integral,1,integral*2),
      'norm_bw':(integral,1,integral*2),
      'norm_cb':(integral,1,integral*2),
      'norm_experf':(integral*0.2,1,integral*2),
      'w_mu':(80,78,82),
      'w_sigma':(8,5,20),
      }

  if not paramset:
    params = dparams['pass'] if isPass else dparams['fail']
  else:
    params = dparams[paramset]

  suffix = '_pass' if isPass else '_fail'
  model = cat.buildModel('prong%i'%nprongs + suffix,mass,model,params)
  if smear:
    smeared_model = cat.smearModel('prong%i_smeared'%nprongs,mass,model,smear,smear_params)
    smeared_model.fitTo(ds,root.RooFit.SumW2Error(False))
  else:
    model.fitTo(ds,root.RooFit.SumW2Error(False))

  Load('Drawers','CanvasDrawer')

  plot = root.CanvasDrawer()
  plot.Ratio(True)
  plot.SetTDRStyle()
  c = root.TCanvas()
  plot.SetCanvas(c)
  plot.SplitCanvas()
  plot.AddCMSLabel()
  pad1 = plot.GetPad1()
  pad1.cd()

  mass.setBins(nbins)
  frame = mass.frame(root.RooFit.Bins(nbins))
  ds.plotOn(frame,root.RooFit.DrawOption('EP'),root.RooFit.MarkerSize(0.8),root.RooFit.Name('MC'))
  model.plotOn(frame,root.RooFit.Name('Fit'))
  frame.GetYaxis().SetTitle('Events/10 GeV')
  frame.SetMaximum(frame.GetMaximum()*1.5)
  frame.Draw()

  pad2 = plot.GetPad2(); pad2.cd()
  frame2 = mass.frame(root.RooFit.Bins(nbins))
  frame2.addObject(frame.pullHist('MC','Fit'),'EP')
  frame2.SetMinimum(-5)
  frame2.SetMaximum(5)
  frame2.GetXaxis().SetTitle('m_{SD} [GeV]')
  frame2.GetYaxis().SetTitle("#frac{Data-Fit}{#sigma_{Data}}");
  frame2.GetYaxis().SetNdivisions(5);
  frame2.GetYaxis().SetTitleSize(15);
  frame2.GetYaxis().SetTitleFont(43);
  frame2.GetYaxis().SetTitleOffset(1.55);
  frame2.GetYaxis().SetLabelFont(43); 
  frame2.GetYaxis().SetLabelSize(15);
  frame2.GetXaxis().SetTitleSize(20);
  frame2.GetXaxis().SetTitleFont(43);
  frame2.GetXaxis().SetTitleOffset(4.);
  frame2.GetXaxis().SetLabelFont(43);
  frame2.GetXaxis().SetLabelSize(15);
  frame2.Draw()

  hZero = root.TH1D('zero','zero',nbins,50,450)
  hZero.SetLineColor(1)
  hZero.Draw('hist same')

  pad1.cd()

  plot.Draw('~/public_html/figs/toptagging/datavalidation/v8/templates/parametric/','prong%i'%nprongs + suffix)

fitMC(3,True ,'CBEXPERF')
fitMC(3,False ,'WGAUSSEXPERF')
fitMC(2,True ,'EXPERF')
fitMC(2,False,'CBEXPERF',paramset='fail2')
fitMC(1,True ,'EXPERF')
fitMC(1,False,'EXPERF')
