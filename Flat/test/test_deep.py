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
from PandaAnalysis.Flat.analysis import deep, deepgen
import PandaAnalysis.Flat.selection as selection
import PandaAnalysis.T3.job_utilities as utils

Load('PandaAnalyzer')

skimmer = root.PandaAnalyzer(debug_level)
a = deepgen()
#a.deepGenGrid = True
#a.deepAntiKtSort = True
#a.deepExC = True
#a.processType=root.kTT
a.dump()
skimmer.SetAnalysis(a)
skimmer.AddPresel(root.GenFatJetSel())
#skimmer.SetPreselectionBit(root.PandaAnalyzer.kGenFatJet)

skimmer.firstEvent=0
skimmer.lastEvent=100
skimmer.isData=False
if skimmer.isData:
    with open(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/certs/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt') as jsonFile:
        payload = json.load(jsonFile)
        for run,lumis in payload.iteritems():
            for l in lumis:
                skimmer.AddGoodLumiRange(int(run),l[0],l[1])
fin = root.TFile.Open(torun)

tree = fin.FindObjectAny("events")
hweights = fin.FindObjectAny("hSumW")
hnpuweights = fin.FindObjectAny("hNPVTrue")
weights = fin.FindObjectAny('weights')
if not hnpuweights:
    hnpuweights = None
if not weights:
    weights = None

skimmer.SetDataDir(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/')
skimmer.Init(tree,hweights,hnpuweights,weights)
skimmer.SetOutputFile(output)

skimmer.Run()
skimmer.Terminate()

#deep_utils.STORE = True
#deep_utils.SAVE = True
#deep_utils.NORM = False
#deep_utils.run_model('testskim_pf_%i.root', 'testskim.root')
