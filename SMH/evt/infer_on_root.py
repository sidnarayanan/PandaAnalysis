#!/usr/bin/env python

from PandaCore.Tools.script import * 
from PandaCore.Tools.root_interface import Selector
import PandaCore.Tools.Functions
from ROOT import gROOT, TFile, TTree, gDirectory, AddressOf, TLorentzVector
import json 
import numpy as np
import os

gROOT.ProcessLine(
"struct TreeStruct {\
Float_t classifier;\
}")

##############################################################################

args = parse('--ifile', '--h5', '--json', '--bname')

from keras.models import load_model 
from ROOT import TreeStruct

with open(args.json) as fjson:
  payload = json.load(fjson)
  features = payload['features']

f = root.TFile(args.ifile)
t = f.Get('events')

treestruct = TreeStruct()
nent = int(t.GetEntries())

#if t.FindBranch(args.bname):
#  t.SetBranchStatus(args.bname, 0)

output = args.ifile + "_tmp"
ofile = root.TFile(output,"RECREATE");
otree = t.CloneTree()

s = Selector()
s.read_tree(t, branches=features)

model = load_model(args.h5)
print 'Loaded model...'

val = array('f', [-1])
o_val = otree.Branch(args.bname, val, args.bname+'/F')

for i in range(t.GetEntriesFast()):
  val[0] = -1
  x = np.array([float(s[f][i]) for f in features])
  x = x.reshape((1,x.shape[0]))
#  print model.predict(x)
  val[0] = model.predict(x)[0,0]
  o_val.Fill()
          

ofile.WriteTObject(otree, 'events')
ofile.Close()
os.system("mv -f %s %s" % (output, args.ifile))

