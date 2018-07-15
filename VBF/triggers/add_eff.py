#!/usr/bin/env python

from sys import argv, exit
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--histfile',type=str)
parser.add_argument('--histname',type=str)
parser.add_argument('--outbranch',type=str,default='sf_l1')
parser.add_argument('outfile',metavar='outfile',type=str,nargs='+',help='input file(s) to process')
args = parser.parse_args()

argv = []

import ROOT as root
from PandaCore.Utils.load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions 

Load('PandaCoreTools')
 
hba = root.H2BranchAdder()
hfile = root.TFile.Open(args.histfile)
ptformula = 'jot%iPt'
etaformula = 'fabs(jot%iEta)'
hbranch = args.outbranch+'_jot%i'
h = hfile.Get(args.histname)
hba.setH(h)

fba = root.FormulaBranchAdder()
fba.newBranchName = args.outbranch
fba.formula = '(1-{0}_jot1)*(1-{0}_jot2)'.format(args.outbranch)
#fba.formula = '((1-sf_l1_jot1) * (1-sf_l1_jot2))'

for fpath in args.outfile:
    f = root.TFile.Open(fpath,'update')
    t = f.Get('events')
    for i in [1,2]:
        hba.formulaX = ptformula%i
        hba.formulaY = etaformula%i 
        hba.newBranchName = hbranch%i 
        hba.addBranch(t)
    f.WriteTObject(t, 'events', 'overwrite')
    f.Close()
    f = root.TFile.Open(fpath,'update')
    t = f.Get('events')
    fba.addBranch(t)
    f.WriteTObject(t, 'events', 'overwrite')
    f.Close()
