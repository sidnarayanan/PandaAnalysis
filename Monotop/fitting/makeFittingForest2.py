#!/usr/bin/env python
from re import sub
from sys import argv,exit
from os import path,getenv
from glob import glob
import argparse
from pprint import pprint

parser = argparse.ArgumentParser(description='make forest')
parser.add_argument('--region',metavar='region',type=str,default=None)
parser.add_argument('--masses',metavar='masses',type=str,default=None) # if masses is set, all couplings for this mass point are run
args = parser.parse_args()
out_region = args.region
region = out_region.split('_')[0]
if region=='test':
    is_test = True 
    region = 'signal'
else:
    is_test = False
masses = args.masses.strip()
if out_region!='signal_vector' and out_region!='signal_scalar':
    masses = None
sname = argv[0]

argv=[]
import PandaAnalysis.Flat.fitting_forest as forest 
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions # kinematics
import PandaAnalysis.Monotop.CombinedBVetoSelection as sel
#import PandaAnalysis.Monotop.CombinedSelection as sel

if masses:
    PDebug(sname,'Looking for point="%s"'%(masses))

basedir = getenv('PANDA_FLATDIR')+'/'
lumi = 35900

def f(x):
    return basedir + x + '.root'

def shift_btags(additional=None):
    shifted_weights = {}
    #if not any([x in region for x in ['signal','top','w']]):
    #    return shifted_weights 
    for shift in ['BUp','BDown','MUp','MDown']:
        for cent in ['sf_btag','sf_sjbtag']:
            shiftedlabel = ''
            if 'sj' in cent:
                shiftedlabel += 'sj'
            if 'B' in shift:
                shiftedlabel += 'btag'
            else:
                shiftedlabel += 'mistag'
            if 'Up' in shift:
                shiftedlabel += 'Up'
            else:
                shiftedlabel += 'Down'
            weight = sel.weights[region+'_'+cent+shift]%lumi
            if additional:
                weight = tTIMES(weight,additional)
            shifted_weights[shiftedlabel] = weight
    return shifted_weights


vmap = {'top_ecf_bdt':'top_ecf_bdt'}
mc_vmap = {'genBosonPt':'genBosonPt'}
if region in ['signal','test']:
    u,uphi, = ('pfmet','pfmetphi')
elif 'photon' in region:
    u,uphi = ('pfUAmag','pfUAphi')
elif 'single' in region:
    u,uphi = ('pfUWmag','pfUWphi')
elif 'di' in region:
    u,uphi = ('pfUZmag','pfUZphi')
vmap['met'] = 'min(%s,999.9999)'%u 
weights = {'nominal' : sel.weights[region]%lumi}

if masses:
    if out_region=='signal_vector':
        with open(getenv('CMSSW_BASE')+'/src/PandaAnalysis/Monotop/fitting/signal_weights_all.dat') as fweights:
            for line in fweights:
                l = line.strip()
                if l!='nominal':
                    vmap[l.replace('rw_','').replace('_nlo','')] = l
    vmap.update(shift_btags())
else:
    weights.update(shift_btags())


# pprint(vmap)
# pprint(weights)
# exit(1)

factory = forest.RegionFactory(name = region if not(is_test) else 'test',
                               cut = sel.cuts[region],
                               variables = vmap, 
                               mc_variables = mc_vmap, 
                               mc_weights = weights)



if is_test:
    factory.add_process(f('Diboson'),'Diboson')
elif region=='photon':
    factory.add_process(f('GJets'),'Pho')
    factory.add_process(f('SinglePhoton'),'Data',is_data=True,extra_cut=sel.triggers['pho'])
    factory.add_process(f('SinglePhoton'),'QCD',is_data=True,
                        extra_weights='sf_phoPurity',extra_cut=sel.triggers['pho'])
#elif out_region not in ['signal_scalar','signal_vector','signal_thq','signal_stdm']:
elif not any([x in out_region for x in ['scalar','vector','thq','stdm']]):
    factory.add_process(f('ZtoNuNu'),'Zvv')
    factory.add_process(f('ZJets'),'Zll')
    factory.add_process(f('WJets'),'Wlv')
    factory.add_process(f('TTbar'),'ttbar')
    factory.add_process(f('SingleTop'),'ST')
    factory.add_process(f('Diboson'),'Diboson')
    factory.add_process(f('QCD'),'QCD')
    if 'electron' in region:
        factory.add_process(f('SingleElectron'),'Data',is_data=True,extra_cut=sel.triggers['ele'])
    else:
        factory.add_process(f('MET'),'Data',is_data=True,extra_cut=sel.triggers['met'])
elif 'vector' in out_region:
    signal_files = glob(basedir+'/Vector*root')
    for f in signal_files:
        fname = f.split('/')[-1].replace('.root','')
        signame = fname
        replacements = {
            'Vector_MonoTop_NLO_Mphi-':'',
            '_gSM-0p25_gDM-1p0_13TeV-madgraph':'',
            '_Mchi-':'_',
        }
        for k,v in replacements.iteritems():
            signame = signame.replace(k,v)
        if masses:
            if fname!=masses:
                continue
            masses = signame
        PDebug(sname,"Opening "+f)
        factory.add_process(f,signame)
elif 'scalar' in out_region:
    signal_files = glob(basedir+'/Scalar*root')
    for f in signal_files:
        fname = f.split('/')[-1].replace('.root','')
        signame = fname
        replacements = {
            'Scalar_MonoTop_LO_Mphi-':'',
            '_13TeV-madgraph':'',
            '_Mchi-':'_',
        }
        for k,v in replacements.iteritems():
            signame = signame.replace(k,v)
        if masses:
            if fname!=masses:
                continue
            masses = signame
        factory.add_process(f,'scalar_'+signame)
elif 'thq' in out_region:
    factory.add_process(f('thq'),'thq')
elif 'stdm' in out_region:
    for m in [300,500,1000]:
        factory.add_process(f('ST_tch_DM-scalar_LO-%i_1-13_TeV'%m),'stdm_%i'%m)
        factory.add_process(f('TT_DM-scalar_LO-%i_1-13_TeV'%m),'ttdm_%i'%m)


if is_test:
    out_region = 'test'
if masses:
    out_region += '_'+masses
    factory.run(basedir+'/fitting/signals_3/fittingForest_%s.root'%out_region)
else:
    factory.run(basedir+'/fitting/fittingForest_%s.root'%out_region)

