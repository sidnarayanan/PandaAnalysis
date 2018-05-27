#!/usr/bin/env python

import argparse
from sys import argv
from os import getenv

VERBOSE=True
parser = argparse.ArgumentParser(description='fit stuff')
parser.add_argument('--outdir',metavar='outdir',type=str,default=None)
parser.add_argument('--disc',metavar='disc',type=str,default='top_ecf_bdt')
parser.add_argument('--cut',metavar='cut',type=float,default=0.1)
parser.add_argument('--useMass',metavar='useMass',type=str,default='FALSE')
args = parser.parse_args()
sname = argv[0]
argv=[]


from math import sqrt
import ROOT as root
from PandaCore.Tools.Misc import *
import PandaAnalysis.Monotop.NoMassPFSelection as sel

basedir = getenv('PANDA_FLATDIR')
fdata = root.TFile(basedir+'/MET.root')
tdata = fdata.Get('events')
fmc = root.TFile(basedir+'/ZJets.root')
tmc = fmc.Get('events')

cut = removeCut(sel.cuts['dimuon'],'fj1MSD')
weight = sel.weights['dimuon']%36600

c = root.TCanvas()

def getIntAndErr(t,cut,weight):
  hone = root.TH1F('h1','h1',1,-2,2)
  hone.Sumw2()
  t.Draw('1>>h1',tTIMES(cut,weight))
  val = hone.GetBinContent(1)
  err = hone.GetBinError(1)
  return val,err

def calcEffAndErr(p,perr,f,ferr):
  eff_ = p/(p+f)
  err_ = pow( perr * f / pow(p+f,2) , 2 )
  err_ += pow( ferr * p / pow(p+f,2) , 2 )
  err_ = sqrt(err_)
  return eff_,err_

def calcSFAndErr(d,derr,m,merr):
  sf_ = d/m
  err_ = pow(derr/m,2)
  err_ += pow(d*merr/pow(m,2),2)
  err_ = sqrt(err_)
  return sf_,err_

if args.useMass.upper()=='FALSE':
  if 'tau' in args.disc:
    passcut='%s<%f'
    failcut='%s>%f'
  else:
    passcut='%s>%f'
    failcut='%s<%f'
else:
  if 'tau' in args.disc:
    passcut='%s<%f && fj1MSD>110 && fj1MSD<210'
    failcut='!(%s<%f && fj1MSD>110 && fj1MSD<210)'
  else:
    passcut='%s>%f && fj1MSD>110 && fj1MSD<210'
    failcut='!(%s>%f && fj1MSD>110 && fj1MSD<210)'

# first get the yields for the tagging cut
datapass,datapasserr = getIntAndErr(tdata,tAND(cut,passcut%(args.disc,args.cut)),'1')
datafail,datafailerr = getIntAndErr(tdata,tAND(cut,failcut%(args.disc,args.cut)),'1')
mcpass,mcpasserr = getIntAndErr(tmc,tAND(cut,passcut%(args.disc,args.cut)),weight)
mcfail,mcfailerr = getIntAndErr(tmc,tAND(cut,failcut%(args.disc,args.cut)),weight)

if VERBOSE:
	logger.debug(sname,'%10s %10s %10s %10s'%('pass','fail','passerr','failerr'))
	logger.debug(sname,'%10f %10f %10f %10f'%(datapass,datafail,datapasserr,datafailerr))
	logger.debug(sname,'%10f %10f %10f %10f'%(mcpass,mcfail,mcpasserr,mcfailerr))

mceff,mcerr = calcEffAndErr(mcpass,mcpasserr,mcfail,mcfailerr)
dataeff,dataerr = calcEffAndErr(datapass,datapasserr,datafail,datafailerr)
tagsf, tagsferr = calcSFAndErr(dataeff,dataerr,mceff,mcerr)

logger.info(sname,'SF(%s) = %.4g +/- %.3g'%(passcut%(args.disc,args.cut),tagsf,tagsferr))
