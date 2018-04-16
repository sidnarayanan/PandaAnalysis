#!/usr/bin/env python
from re import sub
from sys import argv,exit
from os import path,getenv
from glob import glob
import argparse
from pprint import pprint

parser = argparse.ArgumentParser(description='make forest')
parser.add_argument('--region',metavar='region',type=str,default=None)
args = parser.parse_args()
out_region = args.region
region = out_region.split('_')[0]
if region=='test':
    is_test = True 
    region = 'signal'
else:
    is_test = False
sname = argv[0]

argv=[]
import PandaAnalysis.Flat.fitting_forest as forest 
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions # kinematics
import PandaAnalysis.VBF.PandaSelection as sel
import ROOT as root

basedir = getenv('PANDA_FLATDIR')+'/'
lumi = 35800

def f(x):
    return basedir + x + '.root'

# variables to import
vmap = {}
mc_vmap = {'genBosonPt':'genBosonPt'}
if region in ['signal','test']:
    u,uphi, = ('pfmet','pfmetphi')
elif 'photon' in region:
    u,uphi = ('pfUAmag','pfUAphi')
elif 'single' in region:
    u,uphi = ('pfUWmag','pfUWphi')
elif 'di' in region:
    u,uphi = ('pfUZmag','pfUZphi')
vmap['met'] = u 
vmap['mjj'] = 'jot12Mass'
vmap['deta'] = 'jot12DEta'
vmap['dphi'] = 'jot12DPhi'
vmap['eta1'] = 'jot1Eta'
vmap['eta2'] = 'jot2Eta'
weights = {'nominal' : sel.weights[region]%lumi}


# build the factory
factory = forest.RegionFactory(name = region if not(is_test) else 'test',
                               cut = sel.cuts[region],
                               variables = vmap, 
                               mc_variables = mc_vmap, 
                               mc_weights = weights)


# create some TChains for the V+jets
tAllW = root.TChain('events')
for f_ in ['WJets','WJets_EWK']:
    tAllW.AddFile(f(f_))
tAllZvv = root.TChain('events')
for f_ in ['ZtoNuNu','ZtoNuNu_EWK']:
    tAllZvv.AddFile(f(f_))
tAllZll = root.TChain('events')
for f_ in ['ZJets','ZJets_EWK']:
    tAllZll.AddFile(f(f_))

if is_test:
    factory.add_process(f('Diboson'),'Diboson')
else:
    # add the V+jets first
    factory.add_process(f('ZtoNuNu'),'Zvv')
    factory.add_process(f('ZtoNuNu'+'_EWK'),'ewkZvv')
    factory.add_process(tAllZvv,'allZvv')
    factory.add_process(f('ZJets'),'Zll')
    factory.add_process(f('ZJets'+'_EWK'),'ewkZll')
    factory.add_process(tAllZll,'allZll')
    factory.add_process(f('WJets'),'Wlv')
    factory.add_process(f('WJets_EWK'),'ewkWlv')
    factory.add_process(tAllW,'allWlv')
    # other backgrounds
    factory.add_process(f('TTbar'),'ttbar')
    factory.add_process(f('SingleTop'),'ST')
    factory.add_process(f('Diboson'),'Diboson')
    factory.add_process(f('QCD'),'QCD')
    # data
    if 'electron' in region:
        factory.add_process(f('SingleElectron'),'Data',is_data=True,extra_cut=sel.triggers['ele'])
    else:
        factory.add_process(f('MET'),'Data',is_data=True,extra_cut=sel.triggers['met'])
    # signals
    if 'signal' in region:
        factory.add_process(f('ggFHinv_m125'),'GGF_H125')
        for m in [1000, 110, 125, 150, 200, 300, 500, 600, 800 ]:
            factory.add_process(f('vbfHinv_m%i'%m),'VBF_H%i'%m)

factory.run(basedir+'/fitting/fittingForest_%s.root'%out_region)
