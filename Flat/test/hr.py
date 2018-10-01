#!/usr/bin/env python

from re import sub
from sys import argv,exit
from os import system,getenv
import json

debug_level = 0
torun = argv[1]
output = 'testhr.root'
if len(argv)>2:
    debug_level = int(argv[2])
    if len(argv)>3:
        output = argv[3]

argv = []

import ROOT as root
from PandaCore.Utils.load import *
from PandaAnalysis.Flat.analysis import *
import PandaAnalysis.T3.job_utilities as utils
import PandaAnalysis.Tagging.cfg_v8 as tagcfg

Load('PandaAnalysisFlat')

data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'

class BDTAdder(object):
    def __init__(self):
        Load('TMVABranchAdder')
        self.ba = root.TMVABranchAdder()
        self.ba.defaultValue = -1.2
        self.ba.presel = 'clf_ECFN_2_4_20>0'
        for v in tagcfg.variables:
            self.ba.AddVariable(v[0],v[2].replace('fj1','clf_'))
        for v in tagcfg.formulae:
            self.ba.AddFormula(v[0],v[2].replace('fj1','clf_'))
        for s in tagcfg.spectators:
            self.ba.AddSpectator(s[0])
        self.ba.BookMVA('clf_top_ecf_bdt',data_dir+'/trainings/top_ecfbdt_v8_BDT.weights.xml')
    def __call__(self, fname='output.root', tname='events'):
        # now run the BDT
        self.ba.treename = tname
        self.ba.RunFile(fname)

add_bdt = BDTAdder() #backwards compatability

a = analysis("substructure")
a.inpath = torun
a.outpath = 'testhr.root'
a.datapath = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'
a.processType = root.pa.kTop
a.isData = False
utils.set_year(a, 2016)

skimmer = root.pa.HRAnalyzer(a, debug_level)

# skimmer.firstEvent=0
skimmer.lastEvent=10

skimmer.Run()
skimmer.Terminate()

add_bdt(fname='testhr.root')
