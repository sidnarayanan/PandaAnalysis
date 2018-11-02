#!/usr/bin/env python

from sys import argv, exit
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--histfile',type=str)
parser.add_argument('--outbranch',type=str,default='sf_l1')
parser.add_argument('outfile',metavar='outfile',type=str,nargs='+',help='input file(s) to process')
args = parser.parse_args()

argv = []

import ROOT as root
from PandaCore.Utils.load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions 

Load('PandaCoreTools')
 
hba = { 'smu' : root.H2BranchAdder(),
        'jet' : root.H2BranchAdder() }
hfile = root.TFile.Open(args.histfile)
ptformula = 'jotPt[%i]'
etaformula = 'fabs(jotEta[%i])'
hbranch = 'sf_l1%s_jot%i'
hsmu = hfile.Get('h_SingleMuon_spike_finor_pteta')
hjet = hfile.Get('h_JetHT_spike_finor_pteta')
hba['smu'].setH(hsmu)
hba['jet'].setH(hjet)

fba = root.FormulaBranchAdder()
#formula = 'TMath::Max(0,TMath::Min(1,(1-{0}*((sf_l1smu_jot1*(jot1Pt<300)) + (sf_l1jet_jot1*(jot1Pt>300))))*(1-{0}*((sf_l1smu_jot1*(jot1Pt<300)) + (sf_l1jet_jot1*(jot1Pt>300))))))'
formula = '(1-sf_l1smu_jot0)*(1-sf_l1smu_jot1)'
#fba.formula = '(1-{0}_jot1)*(1-{0}_jot2)'.format(args.outbranch)
#fba.formula = '((1-sf_l1_jot1) * (1-sf_l1_jot2))'

for fpath in args.outfile:
    f = root.TFile.Open(fpath,'update')
    t = f.Get('events')
    for i in [0,1]:
        for pd in ['smu','jet']:
            hba[pd].formulaX = ptformula%i
            hba[pd].formulaY = etaformula%i 
            hba[pd].newBranchName = hbranch%(pd,i)
            hba[pd].addBranch(t)
    f.WriteTObject(t, 'events', 'overwrite')
    f.Close()
    f = root.TFile.Open(fpath,'update')
    t = f.Get('events')
    fba.newBranchName = args.outbranch
    fba.formula = formula.format(1)
    fba.addBranch(t)
    '''
    fba.newBranchName = args.outbranch+'_Up'
    fba.formula = formula.format(1.2)
    fba.addBranch(t)
    fba.newBranchName = args.outbranch+'_Down'
    fba.formula = formula.format(0.8)
    fba.addBranch(t)
    '''
    f.WriteTObject(t, 'events', 'overwrite')
    f.Close()
