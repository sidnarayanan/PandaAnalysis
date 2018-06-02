#!/usr/bin/env python

from sys import exit,argv
import argparse

parser = argparse.ArgumentParser(description='skim a file')
parser.add_argument('--input',type=str)
parser.add_argument('--output',type=str,default=None)
parser.add_argument('--cut',type=str)
args = parser.parse_args()

argv=[]

from PandaCore.Utils.load import *
import ROOT as root
Load('Cutter')

c = root.Cutter()
c.treeName = 'events'

if not args.output:
  args.output = args.input.replace('.root','_slim.root')

c.Cut(args.input,args.output,args.cut)
