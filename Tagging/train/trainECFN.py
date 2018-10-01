#!/usr/bin/env python

import ROOT as root
from PandaCore.Tools.Load import *
from PandaCore.Tools.Misc import *
from os import getenv

import cfgInvar as cfg
#import makeInvariantRatios as cfg

Load('Learning','TMVATrainer')

workdir = getenv('SCRAMJETFLAT')

trainer = root.TMVATrainer('top_ecfbdt_50',workdir+'training/')
trainer.treename = 'puppiCA15'
trainer.sigweight = 'ptweight*normalizedWeight'
trainer.bgweight = 'ptweight_analytic*normalizedWeight'
sanitycut='tau32>0 && pt<1000 && mSD>110 && mSD<210'
#trainer.bgcut = tAND(sanitycut,'eventNumber%15==0')
trainer.bgcut =sanitycut
trainer.sigcut = tAND(sanitycut,'matched==1 && gensize<1.2')

for v in cfg.variables:
  trainer.AddVariable(v[0],v[1])

for v in cfg.formulae:
  trainer.AddVariable(v[0],v[1])

for s in cfg.spectators:
  trainer.AddSpectator(s[0],s[1])

trainer.SetFiles(workdir+'/ZpTT.root',workdir+'/QCD_evt10.root')
trainer.BookBDT()
trainer.TrainAll()

