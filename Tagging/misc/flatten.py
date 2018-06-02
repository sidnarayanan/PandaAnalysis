#!/usr/bin/env python

#TODO: replace TTree.Draw with read_tree and fill_hist

from sys import argv
from os import getenv
which = argv[1]
argv = []

import ROOT as root
from PandaCore.Utils.load import *
from PandaCore.Tools.Misc import *

Load('BranchAdder')
 
weight_name = '1' # 'normalizedWeight'

def gethisto(tree,formula,binlo,binhi,additionalcut=None):
  nbins=50
  h = root.TH1D('h','h',nbins,binlo,binhi)
  if additionalcut:
    tree.Draw(formula+'>>h',tTIMES(weight_name,additionalcut))
  else:
    tree.Draw(formula+'>>h',weight_name)

  for iB in xrange(1,nbins+1):
    val = h.GetBinContent(iB);
    if val:
      h.SetBinContent(iB,1./val)
    else:
      print formula,h.GetBinCenter(iB)

  return h


def addbranchesFormula(fpath,additionalcut=None):
  fin = root.TFile(fpath,'UPDATE')
  jets = fin.Get('events')
  ba = root.BranchAdder()

  ba.formula = '1./(250*TMath::Exp((232-fj1Pt)*0.0396)+200*TMath::Exp((235-fj1Pt)*0.0157)+TMath::Exp((583.-fj1Pt)*0.00672))'
  ba.newBranchName = 'ptweight_qcd'
  ba.AddBranchFromFormula(jets)

  fin.WriteTObject(jets,'events','Overwrite')
  fin.Close()

def addbranches(fpath,additionalcut=None):
  fin = root.TFile(fpath,'UPDATE')
  jets = fin.Get('events')
  ba = root.BranchAdder()

  hpt = gethisto(jets,'fj1Pt',250,1000,additionalcut)
  ba.formula = 'fj1Pt'
  ba.newBranchName = 'ptweight_fixed'
  ba.AddBranchFromHistogram(jets,hpt)

  fin.WriteTObject(jets,'events','Overwrite')
  fin.Close()

basedir = getenv('PANDA_FLATDIR')

additionalcut = None
if any([x in which for x in ['Zp','Top']]):
  additionalcut = 'fj1IsMatched==1&&fj1GenSize<1.44'
addbranches(basedir+'/'+which+'.root',additionalcut)
#addbranchesFormula(basedir+'/'+which+'.root')
