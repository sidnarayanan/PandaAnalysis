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

a = breg()
a.bjetBDTReg = False
a.bjetDeepReg = False
a.inpath = torun
a.outpath = 'testskim.root'
a.datapath = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'
a.isData = False
a.processType = root.pa.kTT 
utils.set_year(a, 2017)

skimmer = root.pa.PandaAnalyzer(a, debug_level)

#skimmer.firstEvent=0
skimmer.lastEvent=10
if a.isData:
    utils.add_json(skimmer)

#skimmer.AddPresel(root.pa.LowGenBosonPtSel())
#skimmer.AddPresel(root.pa.VqqHbbSel())
#skimmer.AddPresel(root.pa.TriggerSel())

skimmer.Run()
skimmer.Terminate()
