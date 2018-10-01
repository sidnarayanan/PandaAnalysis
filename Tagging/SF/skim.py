#!/usr/bin/env python

from sys import argv 
from os import getenv,system
me = argv[0]
cat = argv[1]
radius = argv[2]
argv = []

from PandaCore.Utils.load import Load 
import PandaAnalysis.Tagging.TnPSel as sel
from PandaCore.Tools.Misc import tAND
import ROOT as root 

Load('Cutter')

flatdir = getenv('PANDA_FLATDIR')+'/'
outdir = flatdir + 'tnp_'+radius+'/'
system('mkdir -p '+outdir)

fins = {
        'Data' : ['SingleMuon'],
        '1' : ['WJets'],
        '2' : ['Diboson','SingleTop','TTbar'],
        '3' : ['SingleTop','TTbar']
        }
cuts = {
        'Data' : '1==1',
        '1' : '1==1',
        '2' : '(fjIsMatched==0||fjGenSize>%s)'%radius,
        '3' : '(fjIsMatched==1&&fjGenSize<%s)'%radius
       }

c = root.Cutter()
for f in fins[cat]:
    c.Cut(flatdir+f+'.root',
          outdir+f+'.root',
          tAND(cuts[cat], sel.cuts['tag']))

unmerged = [outdir+f+'.root' for f in fins[cat]]
system('hadd -f '+outdir+cat+'.root '+' '.join(unmerged))
[system('rm %s%s.root'%(outdir, f)) for f in fins[cat]]
