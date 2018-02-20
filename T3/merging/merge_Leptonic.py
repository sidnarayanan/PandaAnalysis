#!/usr/bin/env python

from array import array
from glob import glob
from re import sub
from sys import argv,exit
from os import environ,system,path

sname = argv[0]
arguments = [x for x in argv[1:]] # deep copy
argv=[]

import ROOT as root
from PandaCore.Tools.process import *
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import Load

Load('Normalizer')

pds = {}
for k,v in processes.iteritems():
    if v[1]=='MC':
        pds[v[0]] = (k,v[2])  
    else:
        pds[v[0]] = (k,-1)

VERBOSE=True

user = environ['USER']
system('mkdir -p /tmp/%s/split'%user) # tmp dir
system('mkdir -p /tmp/%s/merged'%user) # tmp dir

inbase = environ['SUBMIT_OUTDIR']
outbase = environ['PANDA_FLATDIR']

if VERBOSE:
    suffix = ''
else:
    suffix = ' > /dev/null '

def hadd(inpath,outpath):
    if type(inpath)==type('str'):
        infiles = glob(inpath)
        PInfo(sname,'hadding %s into %s'%(inpath,outpath))
        cmd = 'hadd -k -ff -n 100 -f %s %s %s'%(outpath,inpath,suffix)
        system(cmd)
        return
    else:
        infiles = inpath
    if len(infiles)==0:
        PWarning(sname,'nothing hadded into',outpath)
        return
    elif len(infiles)==1:
        cmd = 'cp %s %s'%(infiles[0],outpath)
    else:
        cmd = 'hadd -k -ff -n 100 -f %s '%outpath
        for f in infiles:
            if path.isfile(f):
                cmd += '%s '%f
    if VERBOSE: PInfo(sname,cmd)
    system(cmd+suffix)

def normalizeFast(fpath,opt):
    xsec=-1
    if type(opt)==type(1.) or type(opt)==type(1):
        xsec = opt
    else:
        try:
            xsec = processes[proc][2]
        except KeyError:
            for k,v in processes.iteritems():
                if proc in k:
                    xsec = v[2]
    if xsec<0:
        PError(sname,'could not find xsec, skipping %s!'%opt)
        return
    xsec = xsec * 1000.0;
    PInfo(sname,'normalizing %s (%s) in fb ...'%(fpath,opt))
    n = root.Normalizer();
    n.NormalizeTree(fpath,xsec)

def merge(shortnames,mergedname):
    for shortname in shortnames:
        if 'monotop' in shortname:
            pd = shortname
            xsec = 1
        elif 'Vector' in shortname:
            tmp_ = shortname
            replacements = {
                'Vector_MonoTop_NLO_Mphi-':'',
                '_gSM-0p25_gDM-1p0_13TeV-madgraph':'',
                '_Mchi-':'_',
                }
            for k,v in replacements.iteritems():
                tmp_ = tmp_.replace(k,v)
            m_V,m_DM = [int(x) for x in tmp_.split('_')]
            params = read_nr_model(m_V,m_DM)
            if params:
                xsec = params.sigma
            else:
                exit(1)
        elif 'Scalar' in shortname:
            tmp_ = shortname
            replacements = {
                'Scalar_MonoTop_LO_Mphi-':'',
                '_13TeV-madgraph':'',
                '_Mchi-':'_',
                }
            for k,v in replacements.iteritems():
                tmp_ = tmp_.replace(k,v)
            m_V,m_DM = [int(x) for x in tmp_.split('_')]
            params = read_r_model(m_V,m_DM)
            if params:
                xsec = params.sigma
            else:
                exit(1)
        elif shortname in pds:
            pd = pds[shortname][0]
            xsec = pds[shortname][1]
        else:
            for shortname_ in [shortname.split('_')[0],shortname.split('_')[-1]]:
                if shortname_ in pds:
                    pd = pds[shortname_][0]
                    xsec = pds[shortname_][1]
                    break
        inpath = inbase+shortname+'_*.root'
        hadd(inpath,'/tmp/%s/split/%s.root'%(user,shortname))
        if xsec>0:
            normalizeFast('/tmp/%s/split/%s.root'%(user,shortname),xsec)
    hadd(['/tmp/%s/split/%s.root'%(user,x) for x in shortnames],'/tmp/%s/merged/%s.root'%(user,mergedname))

d = {
    'test'                          : ['Diboson_ww'],
    'qqZZ'                          : ['ZZTo2L2Nu_13TeV_powheg_pythia8','ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8','ZZTo4L_13TeV_powheg_pythia8'],
    'ggZZ'                          : ['GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8','GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8','GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8','GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8','GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8','GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8','GluGluToContinToZZTo2e2nu_13TeV_MCFM701_pythia8','GluGluToContinToZZTo2mu2nu_13TeV_MCFM701_pythia8',],
    'WZ'                            : ['WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8','WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8','WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8','WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8'],
    'WWSSEWK'                       : ['WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8'],
    'WWSSQCD'                       : ['WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8'],
    'qqWW'                          : ['WWTo2L2Nu_13TeV-powheg'],
    'ggWW'                          : ['GluGluWWTo2L2Nu_MCFM_13TeV'],
    'WWdps'                         : ['WWTo2L2Nu_DoubleScattering_13TeV-pythia8'],
    'VVV'                           : ['WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8','WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8','WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8','ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8','WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8'],
    'TTV'                           : ['TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8','TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8','TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8','TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8','TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8','tZq_ll_4f_13TeV-amcatnlo-pythia8','tZq_nunu_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1'],
    'TT2L'                          : ['TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8'],
    'TW'                            : ['ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1','ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'],
    'WGstar'                        : ['WGstarToLNuEE_012Jets_13TeV-madgraph','WGstarToLNuMuMu_012Jets_13TeV-madgraph'],
    'VG'                            : ['WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8','ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'],
    'H125'                          : ['GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8','VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8','GluGluHToWWTo2L2Nu_M125_13TeV_powheg_JHUgen_pythia8','VBFHToWWTo2L2Nu_M125_13TeV_powheg_JHUgenv628_pythia8','GluGluHToTauTau_M125_13TeV_powheg_pythia8','VBFHToTauTau_M125_13TeV_powheg_pythia8','ttHToNonbb_M125_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8','VHToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8'],
    'DYNJetsToLL'                   : ['DY1JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8','DY2JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8','DY3JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8','DY4JetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8','DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8','DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8','DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8','DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
    'DYJetsToLL_Pt0To50'            : ['DYJetsToLL_Zpt-0To50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'],
    'DYJetsToLL_Pt50To100'          : ['DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'],
    'DYJetsToLL_Pt100To250'         : ['DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'],
    'DYJetsToLL_Pt250To400'         : ['DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'],
    'DYJetsToLL_Pt400To650'         : ['DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'],
    'DYJetsToLL_Pt650ToInf'         : ['DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'],
    'DYJetsToTauTau'                : ['DYJetsToTauTau_ForcedMuEleDecay_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'],
    'DYJetsToLL_M-10to50'           : ['DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
    'DYJetsToLL_M-50_NLO'           : ['DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'],
    'DYJetsToLL_M-50_LO_Pt000To050' : ['DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
    'DYJetsToLL_M-50_LO_Pt100to200' : ['DYJetsToLL_Zpt-100to200_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
    'DYJetsToLL_M-50_LO_Pt200toInf' : ['DYJetsToLL_Zpt-200toInf_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'],
    'DYJetsToLL_POWHEG'             : ['ZToMuMu_NNPDF30_13TeV-powheg_M_50_120','ZToMuMu_NNPDF30_13TeV-powheg_M_120_200','ZToEE_NNPDF30_13TeV-powheg_M_50_120','ZToEE_NNPDF30_13TeV-powheg_M_120_200'],
    'DYJetsToMM_POWHEG'             : ['ZToMuMu_NNPDF30_13TeV-powheg_M_50_120'],
    'DYJetsToEE_POWHEG'             : ['ZToEE_NNPDF30_13TeV-powheg_M_50_120'],
    'WJets'                         : ['WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'],
    'TT'                            : ['TT_TuneCUETP8M2T4_13TeV-powheg-pythia8'],
    'WpWp_EWK'                      : ['WpWpJJ_EWK_TuneCUETP8M1_13TeV-madgraph-pythia8'],
    'WpWp_QCD'                      : ['WpWpJJ_QCD_TuneCUETP8M1_13TeV-madgraph-pythia8'],
    'MET'                           : ['MET'],
    'SingleMuon'                    : ['SingleMuon'],
    'SingleElectron'                : ['SingleElectron'],
    'MuonEG'                        : ['MuonEG'],
    'DoubleEG'                      : ['DoubleEG'],
    'DoubleMuon'                    : ['DoubleMuon'],
    'data_overlaps'                 : ['MuonEG','DoubleMuon','SingleMuon','DoubleEG','SingleElectron'],
    'SinglePhoton'                  : ['SinglePhoton'],
    'QCD'                           : ['QCD_ht100to200','QCD_ht200to300','QCD_ht300to500','QCD_ht500to700','QCD_ht700to1000','QCD_ht1000to1500','QCD_ht1500to2000','QCD_ht2000toinf']
}

args = {}

for pd in arguments:
    if pd in d:
        args[pd] = d[pd]
    else:
        args[pd] = [pd]

for pd in args:
    merge(args[pd],pd)
    system('cp -r /tmp/%s/merged/%s.root %s'%(user,pd,outbase))
    PInfo(sname,'finished with '+pd)

