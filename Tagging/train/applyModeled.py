#!/usr/bin/env python

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
from os import getenv
from sys import argv

import cfgModeled as cfg

Load('Learning','TMVABranchAdder')

workdir = getenv('SCRAMJETFLAT')

ba = root.TMVABranchAdder()
ba.treename='puppiCA15'

for v in cfg.variables:
  ba.AddVariable(v[0],v[0])

for v in cfg.formulae:
  ba.AddFormula(v[0],v[0])

for s in cfg.spectators:
  ba.AddSpectator(s[0])

#ba.BookMVA('top_ecfv7_bdt',workdir+'/training/top_ecfbdt_v7_BDT.weights.xml')
ba.BookMVA('top_ecfv6_bdt',workdir+'/training/top_ecfbdt_v6_BDT.weights.xml')

ba.RunFile(workdir+'/'+argv[1]+'.root')

