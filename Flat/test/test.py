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
from PandaCore.Utils.load import *
from PandaAnalysis.Flat.analysis import *
import PandaAnalysis.T3.job_utilities as utils

Load('PandaAnalyzer')

<<<<<<< HEAD
skimmer = root.PandaAnalyzer(debug_level)
#gghbb = gghbb()
#gghbb.reclusterGen = False
#gghbb.bjetRegression = False
#gghbb.btagSFs = False
#gghbb.deep = True
#gghbb.dump()
a = vbf()
a.processType = root.kZ
skimmer.SetAnalysis(a)

skimmer.firstEvent=0
skimmer.lastEvent=1000
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
=======
a = monotop(True)
a.recalcECF = True
a.varyJESTotal = True

a.inpath = torun
a.outpath = 'testskim.root'
a.datapath = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'

a.isData = False
utils.set_year(a, 2016)
#utils.set_year(a, 2017)

skimmer = root.pa.PandaAnalyzer(a, debug_level)

#skimmer.firstEvent=0
skimmer.lastEvent=100
if a.isData:
    utils.add_json(skimmer)

#skimmer.AddPresel(root.pa.LowGenBosonPtSel())
#skimmer.AddPresel(root.pa.VHbbSel())
#skimmer.AddPresel(root.pa.TriggerSel())
>>>>>>> sid/master

skimmer.Run()
skimmer.Terminate()
