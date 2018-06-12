#!/usr/bin/env python

from sys import argv
from os import getenv
import argparse
parser = argparse.ArgumentParser(description='reweight CSV shape for fakes')
parser.add_argument('--reweight',metavar='reweight',type=str)
parser.add_argument('--infile',metavar='infile',type=str)
args = parser.parse_args()
argv=[]

import ROOT as root
from PandaCore.Utils.load import *
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions
Load('Tools','BranchAdder')

freweight = root.TFile.Open(args.reweight)
hreweight = freweight.Get('hratio')

ba = root.BranchAdder()
ba.formula = 'isojet1CSV'
ba.newBranchName = 'tagCSVWeightv2'
ba.cut = 'isojet1CSV>0.5 && isojet1Flav!=0' 

fin = root.TFile(args.infile,'UPDATE')
tin = fin.Get('events')
ba.AddBranchFromHistogram(tin,hreweight)
fin.WriteTObject(tin,'events','Overwrite')
fin.Close()
