#!/usr/bin/env python

from PandaCore.Tools.script import * 
import cfg_v8 as cfg 

Load('PandaCoreLearning')

trainer = root.TMVATrainer('top_ecfbdt_v0',flatdir+'training/')
trainer.treename = 'events'
trainer.sigweight = 'weight'
trainer.bgweight = 'weight'
sanitycut='clf_IsMatched==1 && gen_size<1.2 && gen_pt<1200 && clf_MSD>110 && clf_MSD<210'
trainer.bgcut=sanitycut

for v in cfg.variables:
  trainer.AddVariable(v[0],v[1])

for v in cfg.formulae:
  trainer.AddVariable(v[0],v[1])

for s in cfg.spectators:
  trainer.AddSpectator(s[0],s[1])

trainer.SetFiles(flatdir+'/ZpTT.root',flatdir+'/QCD.root')
#trainer.BookBDT(root.TMVATrainer.kAda)
trainer.BookBDT(root.TMVATrainer.kGradWide)
trainer.TrainAll()

