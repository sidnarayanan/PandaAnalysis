#!/usr/bin/env python

from os import getenv
from sys import argv

#import cfgModeled as cfg
import cfgAll as cfg

sample = argv[1]
argv=[]

import ROOT as root
from PandaCore.Utils.load import *
from PandaCore.Tools.Misc import *

Load('TMVABranchAdder')

workdir = getenv('PANDA_FLATDIR')
trainingdir = '/home/snarayan/home000/store/scramjet/v7/training/'

ba = root.TMVABranchAdder()
ba.treename='events'
ba.defaultValue=-1.2
ba.presel='fj1ECFN_2_4_20>0'

for v in cfg.variables:
  ba.AddVariable(v[0],v[2])

for v in cfg.formulae:
  ba.AddFormula(v[0],v[2])

for s in cfg.spectators:
  ba.AddSpectator(s[0])

#ba.BookMVA('top_ecfv6_bdt',workdir+'/training/top_ecfbdt_v6_BDT.weights.xml')
#ba.BookMVA('top_ecfv7_bdt',workdir+'/training/top_ecfbdt_v7_BDT.weights.xml')
#ba.BookMVA('top_ecfv8_bdt',workdir+'/training/top_ecfbdt_v8_BDT.weights.xml')
ba.BookMVA('top_ecfAll_bdt',trainingdir+'scan50/top_ecfbdt_50_BDT.weights.xml')

ba.RunFile(workdir+'/'+sample+'.root')

