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
from PandaAnalysis.Flat.analysis import vv
import PandaAnalysis.T3.job_utilities as utils

Load('PandaAnalyzer')

analysis = vv(True)
analysis.inpath = torun
analysis.outpath = output
analysis.datapath = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'
analysis.isData = False
utils.set_year(analysis, 2017)
analysis.processType = utils.classify_sample(torun, analysis.isData)

print "Process: ",analysis.processType

skimmer = root.pa.PandaAnalyzer(analysis, debug_level)
skimmer.AddPresel(root.pa.LeptonSel())
skimmer.AddPresel(root.pa.TriggerSel())

skimmer.firstEvent=0
skimmer.lastEvent=-1
skimmer.isData=False
if skimmer.isData:
    with open(getenv('CMSSW_BASE')+'/src/PandaAnalysis/data/certs/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt') as jsonFile:
        payload = json.load(jsonFile)
        for run,lumis in payload.iteritems():
            for l in lumis:
                skimmer.AddGoodLumiRange(int(run),l[0],l[1])
fin = root.TFile.Open(torun)

tree = fin.FindObjectAny("events")
hweights = fin.FindObjectAny("hSumW")
weights = fin.FindObjectAny('weights')
if not weights:
    weights = None

#skimmer.Init(tree,hweights,weights)
#skimmer.SetOutputFile(output)

skimmer.Run()
skimmer.Terminate()
