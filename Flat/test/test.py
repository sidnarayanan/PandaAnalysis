#!/usr/bin/env python

from re import sub
from sys import argv,exit
from os import system,getenv
import json

debug_level = 0
torun = argv[1]
output = 'testskim.root'
if len(argv)>2:
    debug_level = int(argv[2])
    if len(argv)>3:
        output = argv[3]

argv = []

import ROOT as root
from PandaCore.Tools.Load import *
from PandaAnalysis.Flat.analysis import *
import PandaAnalysis.T3.job_utilities as utils

Load('PandaAnalyzer')

skimmer = root.PandaAnalyzer(debug_level)
#gghbb = gghbb()
#gghbb.reclusterGen = False
#gghbb.bjetRegression = False
#gghbb.btagSFs = False
#gghbb.deep = True
#gghbb.dump()
a = monotop()
utils.set_year(a, 2017)
skimmer.SetAnalysis(a)

skimmer.firstEvent=0
skimmer.lastEvent=10
skimmer.isData=True
utils.add_json(skimmer)
fin = root.TFile.Open(torun)

tree = fin.FindObjectAny("events")
hweights = fin.FindObjectAny("hSumW")
weights = fin.FindObjectAny('weights')
if not weights:
    weights = None

skimmer.SetDataDir(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/')
skimmer.Init(tree,hweights,weights)
skimmer.SetOutputFile(output)

skimmer.Run()
skimmer.Terminate()
