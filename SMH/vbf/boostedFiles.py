#!/usr/bin/env python

from sys import exit,argv,stdout
import argparse
from glob import glob
from os import getenv, system
from multiprocessing.pool import ThreadPool

parser = argparse.ArgumentParser(description='skim a file')
parser.add_argument('--cut',type=str,default='fjPt>400')
args = parser.parse_args()

argv=[]

from PandaCore.Tools.Load import *
import ROOT as root
Load('Cutter')

fdir = getenv('PANDA_FLATDIR')
system('mkdir -p %s/boosted/'%fdir)

def cut(f):
  c = root.Cutter()
  c.treeName = 'events'
  fbase = f.split('/')[-1]
  stdout.write(fbase+'\n'); stdout.flush()
  c.Cut(f,fdir+'/boosted/'+fbase,args.cut)


pool = ThreadPool(5)
pool.map(cut, glob(fdir+'/*root'))
