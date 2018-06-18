#!/usr/bin/env python

from re import sub
from sys import argv,exit
from os import system,getenv
import json

debug_level = 0
torun = argv[1]
output = 'testl1.root'
if len(argv)>2:
    debug_level = int(argv[2])
    if len(argv)>3:
        output = argv[3]

argv = []

import ROOT as root
from PandaCore.Utils.load import *
from PandaAnalysis.Flat.analysis import *
import PandaAnalysis.T3.job_utilities as utils

Load('PandaAnalysisFlat')

a = analysis("l1")
a.inpath = torun
a.outpath = output 
a.datapath = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'

skimmer = root.pa.L1Analyzer(a, debug_level)

#skimmer.firstEvent=0
#skimmer.lastEvent=10

skimmer.Run()
skimmer.Terminate()
