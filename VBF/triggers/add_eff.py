#!/usr/bin/env python

from sys import argv
from os import getenv
which = argv[1]
argv = []

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *

Load('BranchAdder')

cuts = {
    'BB' : 'fabs(jot1Eta)<3 && fabs(jot2Eta)<3',
#    'FF' : 'fabs(jot1Eta)>3 && fabs(jot2Eta)>3',
    'BF' : '!(fabs(jot1Eta)<3 && fabs(jot2Eta)<3) && !(fabs(jot1Eta)>3 && fabs(jot2Eta)>3)'
}

nmu = 2

f = root.TFile(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/vbf16/trig/param_nmu%i.root'%nmu)

def add_individual(key):
  fin = root.TFile(which, 'UPDATE')
  t = fin.Get('events')
  ba = root.BranchAdder()
  ba.defaultValue = 1
  ba.formula = 'jot12Mass'
  ba.cut = cuts[key]
  ba.newBranchName = 'sf_metTrig_nmu%i_%s'%(nmu, key)
  ba.AddBranchFromHistogram(t, f.Get('h_jot12Mass_%s'%key))
  fin.WriteTObject(t, 'events', 'Overwrite')
  fin.Close()
  

def combine():
  fin = root.TFile(which, 'UPDATE')
  t = fin.Get('events')
  ba = root.BranchAdder()
  ba.formula = 'sf_metTrig_nmu%i_BB * sf_metTrig_nmu%i_BF'%(nmu, nmu)
  ba.newBranchName = 'sf_metTrig_nmu%i'%nmu
  ba.AddBranchFromFormula(t)
  fin.WriteTObject(t, 'events', 'Overwrite')
  fin.Close()

for k in cuts:
  add_individual(k)
combine()
