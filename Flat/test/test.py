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

a = wlnhbb()
a.bjetDeepReg = True
a.bjetBDTReg = False
#a = breg()
# a.bjetRegTrain = True
a.inpath = torun
a.outpath = 'testskim.root'
a.datapath = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'
a.processType = root.pa.kZ
a.isData = True
utils.set_year(a, 2017)

skimmer = root.pa.PandaAnalyzer(a, debug_level)

#skimmer.firstEvent=0
skimmer.lastEvent=1000
if a.isData:
    utils.add_json(skimmer)

#skimmer.AddPresel(root.pa.LowGenBosonPtSel())
skimmer.AddPresel(root.pa.VHbbSel())
skimmer.AddPresel(root.pa.TriggerSel())

skimmer.Run()
skimmer.Terminate()
