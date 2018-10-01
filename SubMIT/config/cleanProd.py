#!/usr/bin/env python

from glob import glob
from os import stat,getenv
from re import sub
import sys
import argparse
import subprocess
from PandaCore.Tools.Misc import *

def du(path):
  return float(subprocess.check_output(['du','-s', path]).split()[0].decode('utf-8'))

parser = argparse.ArgumentParser(description='clean up duplicates in a single prod campaign')
args = parser.parse_args()

eosproddirs = getenv('PANDA_PROD').split(':')

datasets = []
for d in eosproddirs:
  datasets += glob(d+'/*/*')

toRemove = []
for dataset in datasets:
  versions = [x.split('/')[-1] for x in glob(dataset+'/*')]
  best = '0'
  bestvol = 0
  for v in versions:
    vvol = du(dataset+'/'+v)
    if v > best and vvol>=bestvol:
      # this is a more recent version and is larger
      best = v
      bestvol = vvol
    else:
      logger.warning('cleanProd.py','dataset=%s'%dataset)
      logger.warning('cleanProd.py','\t %s (%.2f) vs %s (%.2f)'%(v,vvol,best,bestvol))
  for v in versions:
    if v!=best:
      toRemove.append(dataset+'/'+v)


for r in toRemove:
  print r
