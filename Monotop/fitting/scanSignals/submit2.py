#!/usr/bin/env python

from os import system,environ,readlink,path
from sys import exit,stdout,argv
from glob import glob
from PandaCore.Tools.condor import classad,htcondor

user = environ['USER']
fittingdir = environ['PANDA_FITTING']
scansdir = fittingdir+'/scans2/'
scramdir = environ['PANDA_FIT']
flatdir  = environ['PANDA_FLATDIR']

argv = []
import ROOT as root

sigfiles = []
nper = 5
force = True 
#nper = 1

# sigfiles += glob(flatdir+'/ST*root')
# coupling_scan = False

sigfiles += glob(flatdir+'/Scalar*root')
coupling_scan = False

#sigfiles += glob(flatdir+'/Vector*root')
#coupling_scan = True

if coupling_scan:
    fweights = open('../signal_weights.dat')
    weights = [x.strip().replace('_nlo','').replace('rw_','') for x in fweights]
else:
    weights = ['nominal']


cfgheader = '%s %s'%(fittingdir,scansdir)
cfglines = []
mask = {}

for ff in sigfiles:
    f = ff.split('/')[-1].replace('.root','')
    model = ''
    mParams = ''
    mV = ''
    mChi = ''
    couplings = ['']
    if 'Vector' in f:
        replacements = {
            'Vector_MonoTop_NLO_Mphi-':'',
            '_gSM-0p25_gDM-1p0_13TeV-madgraph':'',
            '_Mchi-':'_',
            }
        for k,v in replacements.iteritems():
            f = f.replace(k,v)
        mV,mChi = map(int,f.split('_'))
        model = 'vector'
        for w in weights:
            cfg = '%s %i_%i %s'%(model,mV,mChi,w) 
            if force:
                mask[cfg] = True
            else:
                if len(glob(scansdir+'/'+model+'/'+w+'/higgsCombine*%i_%i*root'%(mV,mChi)))==0:
                    mask[cfg] = True
                else:
                    checkfilepath = glob(scansdir+'/'+model+'/'+w+'/higgsCombine*%i_%i*root'%(mV,mChi))[0]
                    f = root.TFile(checkfilepath)
                    t = f.Get('limit')
                    if not(t) or t.GetEntriesFast()<6:
                        print 'corrupt:',checkfilepath
                        mask[cfg] = True
                    else:
                        mask[cfg] = False
            cfglines.append(cfg)
    else:
        replacements = {
               'Scalar_MonoTop_LO_Mphi-':'',
               '_13TeV-madgraph':'',
               '_Mchi-':'_',
            }
        for k,v in replacements.iteritems():
            f = f.replace(k,v)
        mV,mChi = map(int,f.split('_'))
        model = 'scalar'
        for w in weights:
            cfg = '%s %i_%i %s'%(model,mV,mChi,w) 
            if force:
                mask[cfg] = True
            else:
                if len(glob(scansdir+'/'+model+'/'+w+'/higgsCombine*%i_%i*root'%(mV,mChi)))==0:
                    mask[cfg] = True
                else:
                    checkfilepath = glob(scansdir+'/'+model+'/'+w+'/higgsCombine*%i_%i*root'%(mV,mChi))[0]
                    f = root.TFile(checkfilepath)
                    t = f.Get('limit')
                    if not(t) or t.GetEntriesFast()<6:
                        print 'corrupt:',checkfilepath
                        mask[cfg] = True
                    else:
                        mask[cfg] = False
            cfglines.append(cfg)

to_submit = []
for i in xrange(len(cfglines)):
    cfg = cfglines[i]
    if mask[cfg]:
        to_submit.append(cfg)


# fout = open('submit.cfg','w')
# fout.write(cfgheader+'\n')
# for l in cfglines:
#     fout.write(l+'\n')
# fout.close()

args_template = '%s '%scramdir
args_template += '--template newewk_tmpl.txt '
args_template += '--indir %s '%fittingdir
args_template += '--outdir %s '%scansdir
args_template += '--cfg '

for cfg in to_submit:
    args = args_template
    args += cfg
    print args
