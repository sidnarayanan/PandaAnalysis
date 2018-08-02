#!/usr/bin/env python

from PandaCore.Tools.script import * 
from PandaCore.Tools.root_interface import Selector
import numpy as np 
from array import array 

Load('PandaCoreTools')
args = parse('--infile') 

bins = array('f', np.linspace(250, 1200, 50))
h = root.TH1D('h','h', len(bins)-1, bins)
s = Selector()
f = root.TFile.Open(args.infile, 'update')
t = f.Get('events')
s.read_tree(t, branches=['gen_pt'], cut='clf_IsMatched==1')
h = s.draw(['gen_pt'], hbase=h)

norm = h.Integral() / (len(bins) - 1)
for ib in xrange(1, h.GetNbinsX()+1):
  val = h.GetBinContent(ib)
  if val > 0:
    h.SetBinContent(ib, norm / val)

ba = root.H1BranchAdder()
ba.setH(h)
ba.formulaX = 'gen_pt'
ba.newBranchName = 'weight' 

ba.addBranch(t)

f.WriteTObject(t, 'events', 'overwrite')
f.Close()
