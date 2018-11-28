#!/usr/bin/env python

#TODO: replace TTree.Draw with read_tree and fill_hist
from PandaCore.Tools.script import *
from PandaCore.Tools.root_interface import Selector 

Load('BranchAdder')
 

def gethisto(tree,formula,binlo,binhi,additionalcut=None):
  nbins=100
  h = root.TH1D('h','h',nbins,binlo,binhi)
  s = Selector()
  s.read_tree(tree, branches=[formula], cut=additionalcut)
  h = s.draw(formula, hbase=h)
  h.Smooth()
  for ib in xrange(1, h.GetNbinsX()+1):
    v = h.GetBinContent(ib)
    if v == 0:
      h.SetBinContent(ib,1)
    else:
      h.SetBinContent(ib,1./v)
  return h


def addbranches(fpath,additionalcut=None):
  fin = root.TFile(fpath,'UPDATE')
  jets = fin.Get('events')
  ba = root.H1BranchAdder()

  hpt = gethisto(jets,'gen_pt',0,1000,additionalcut)
  ba.formulaX = 'gen_pt'
  ba.newBranchName = 'ptweight'
  ba.setH(hpt)
  ba.addBranch(jets)

  fin.WriteTObject(jets,'events','Overwrite')
  fin.Close()

basedir = getenv('PANDA_FLATDIR')

addbranches(basedir+'/'+sys.argv[1]+'.root', additionalcut='fabs(gen_eta)<2.5')
