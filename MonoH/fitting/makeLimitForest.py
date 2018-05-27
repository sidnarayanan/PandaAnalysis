#!/usr/bin/env python
from re import sub
from sys import argv,exit
from os import path,getenv
import argparse
parser = argparse.ArgumentParser(description='make forest')
parser.add_argument('--region',metavar='region',type=str,default=None)
toProcess = parser.parse_args().region

argv=[]
import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import *
import PandaCore.Tools.Functions 
import PandaAnalysis.MonoH.Selection as sel
#import PandaAnalysis.MonoH.Selection_mw as sel
#import PandaAnalysis.MonoH.Selection_t21 as sel

Load('PandaAnalysisFlat','LimitTreeBuilder')


subfolder = ''
if toProcess=='signal':
  subfolder = 'preselection/'
elif toProcess=='wmn':
  subfolder = 'cr_w_mu/'
elif toProcess=='wen':
  subfolder = 'cr_w_el/'
elif toProcess=='tm':
  subfolder = 'cr_ttbar_mu/'
elif toProcess=='te':
  subfolder = 'cr_ttbar_el/'
elif toProcess=='zmm':
  subfolder = 'cr_dimuon/'
elif toProcess=='zee':
  subfolder = 'cr_dielectron/'
elif toProcess=='pho':
  subfolder = 'cr_gamma/'
   

#forestDir = '/scratch3/bmaier/flat/'
#baseDir = '/scratch3/bmaier/flat/' + subfolder
forestDir = '/local/bmaier/flat/'
baseDir = '/local/bmaier/flat/' + subfolder
lumi = 12918.

factory = root.LimitTreeBuilder()
if toProcess:
  factory.SetOutFile(forestDir+'/limits/limitForest_%s.root'%toProcess)
else:
  factory.SetOutFile(forestDir+'/limits/limitForest_all.root')

def dataCut(basecut,trigger):
  # return tAND('metFilter==1',tAND(trigger,basecut))
  return tAND(trigger,basecut)

def getTree(fpath):
  fIn = root.TFile(baseDir+fpath+'.root')
  tIn = fIn.Get('events')
  return tIn,fIn

def enable(regionName):
  if toProcess:
    return (toProcess==regionName)
  else:
    return True

def shiftBtags(label,tree,varmap,cut,basecut):
  # basecut is really the baseweight 
  ps = []
  for shift in ['BUp','BDown','MUp','MDown']:
    for cent in ['sf_btag','sf_sjbtag']:
      shiftedlabel = '_'
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
      weight = sel.weights[basecut+'_'+cent+shift]%lumi
      logger.info('makeLimitForest.py',weight)
      shiftedProcess = root.Process(label,tree,varmap,cut,weight)
      shiftedProcess.syst = shiftedlabel
      ps.append(shiftedProcess)
  return ps


# input
tZll,fZll = getTree('ZJets')
tZvv,fZvv = getTree('ZtoNuNu')
tWlv,fWlv = getTree('WJets')
tPho,fPho = getTree('GJets')
tTTbar,fTT = getTree('TTbar')
tVV,fVV = getTree('Diboson')
tQCD,fQCD = getTree('QCD')
tST,fST = getTree('SingleTop')
tMET,fMET = getTree('MET')
tZH,fZH = getTree('ZHbb')
tSingleEle,fSEle = getTree('SingleElectron')
tSinglePho,fSPho = getTree('SinglePhoton')
tZpA0_1000,fZpA0_1000 = getTree('ZpA0-1000-300')
#tZpA0_1200,fZpA0_1200 = getTree('ZpA0-1200-300')
tZpA0_1400,fZpA0_1400 = getTree('ZpA0-1400-300')
tZpA0_1700,fZpA0_1700 = getTree('ZpA0-1700-300')
#tZpA0_2000,fZpA0_2000 = getTree('ZpA0-2000-300')
tZpA0_2500,fZpA0_2500 = getTree('ZpA0-2500-300')


factory.cd()
regions = {}
processes = {}

vm = root.VariableMap()
vm.AddVar('met','puppimet')
vm.AddVar('genBosonPt','genBosonPt')
vmW = root.VariableMap()
vmW.AddVar('met','UWmag')
vmW.AddVar('genBosonPt','genBosonPt')
vmZ = root.VariableMap()
vmZ.AddVar('met','UZmag')
vmZ.AddVar('genBosonPt','genBosonPt')
vmA = root.VariableMap()
vmA.AddVar('met','UAmag')
vmA.AddVar('genBosonPt','genBosonPt')


# test region
if enable('test'):
  regions['test'] = root.Region('test')
  cut = sel.cuts['signal']
  weight = '%f*%s'%(lumi,sel.weights['signal'])
  processes['test'] = [
    root.Process('Data',tMET,vm,dataCut(cut,sel.triggers['met']),'1'),
    root.Process('Diboson',tVV,vm,cut,weight),
  ]

  for p in processes['test']:
    regions['test'].AddProcess(p)
  factory.AddRegion(regions['test'])

# signal region
if enable('signal'):
  regions['signal'] = root.Region('signal')
  cut = sel.cuts['signal']
  weight = '%s'%(sel.weights['signal']%lumi)
  logger.info('makeLimitForest.py',cut)
  logger.info('makeLimitForest.py',weight)
  processes['signal'] = [
    root.Process('Data',tMET,vm,dataCut(cut,sel.metTrigger),'1'),
    root.Process('Zvv',tZvv,vm,cut,weight),
    root.Process('Wlv',tWlv,vm,cut,weight),
    root.Process('Zll',tZll,vm,cut,weight),
    root.Process('ttbar',tTTbar,vm,cut,weight),
    root.Process('ST',tST,vm,cut,weight),
    root.Process('Diboson',tVV,vm,cut,weight),
    root.Process('QCD',tQCD,vm,cut,weight),
    root.Process('ZH',tZH,vm,cut,weight),
    root.Process('ZpA0_1000',tZpA0_1000,vm,cut,weight),
    #root.Process('ZpA0_1200',tZpA0_1200,vm,cut,weight),
    root.Process('ZpA0_1400',tZpA0_1400,vm,cut,weight),
    root.Process('ZpA0_1700',tZpA0_1700,vm,cut,weight),
    #root.Process('ZpA0_2000',tZpA0_2000,vm,cut,weight),
    root.Process('ZpA0_2500',tZpA0_2500,vm,cut,weight),

  ]
  processes['signal'] += shiftBtags('ttbar',tTTbar,vm,cut,'signal')
  processes['signal'] += shiftBtags('Zvv',tZvv,vm,cut,'signal')
  processes['signal'] += shiftBtags('Wlv',tWlv,vm,cut,'signal')
  processes['signal'] += shiftBtags('ST',tST,vm,cut,'signal')
  processes['signal'] += shiftBtags('Zll',tZll,vm,cut,'signal')
  processes['signal'] += shiftBtags('Diboson',tVV,vm,cut,'signal')
  processes['signal'] += shiftBtags('QCD',tQCD,vm,cut,'signal')
  processes['signal'] += shiftBtags('ZH',tZH,vm,cut,'signal')
  processes['signal'] += shiftBtags('ZpA0_1000',tZpA0_1000,vm,cut,'signal')
#  processes['signal'] += shiftBtags('ZpA0_1200',tZpA0_1200,vm,cut,'signal')
  processes['signal'] += shiftBtags('ZpA0_1400',tZpA0_1400,vm,cut,'signal')
  processes['signal'] += shiftBtags('ZpA0_1700',tZpA0_1700,vm,cut,'signal')
#  processes['signal'] += shiftBtags('ZpA0_2000',tZpA0_2000,vm,cut,'signal')
  processes['signal'] += shiftBtags('ZpA0_2500',tZpA0_2500,vm,cut,'signal')
  for p in processes['signal']:
    regions['signal'].AddProcess(p)
  factory.AddRegion(regions['signal'])

# wmn
if enable('wmn'):
  regions['wmn'] = root.Region('wmn')
  cut = sel.cuts['wmn']
  weight = '%s'%(sel.weights['wmn']%lumi)
  logger.info('makeLimitForest.py',cut)
  logger.info('makeLimitForest.py',weight)
  processes['wmn'] = [
    root.Process('Data',tMET,vmW,dataCut(cut,sel.metTrigger),'1'),
    root.Process('Wlv',tWlv,vmW,cut,weight),
    root.Process('Zll',tZll,vmW,cut,weight),
    root.Process('ttbar',tTTbar,vmW,cut,weight),
    root.Process('ST',tST,vmW,cut,weight),
    root.Process('Diboson',tVV,vmW,cut,weight),
    root.Process('QCD',tQCD,vmW,cut,weight),
  ]
  processes['wmn'] += shiftBtags('ttbar',tTTbar,vmW,cut,'wmn')
  processes['wmn'] += shiftBtags('Wlv',tWlv,vmW,cut,'wmn')
  processes['wmn'] += shiftBtags('ST',tST,vmW,cut,'wmn')
  processes['wmn'] += shiftBtags('Zll',tZll,vmW,cut,'wmn')
  processes['wmn'] += shiftBtags('Diboson',tVV,vmW,cut,'wmn')
  processes['wmn'] += shiftBtags('QCD',tQCD,vmW,cut,'wmn')
  for p in processes['wmn']:
    regions['wmn'].AddProcess(p)
  factory.AddRegion(regions['wmn'])

# wen
if enable('wen'):
  regions['wen'] = root.Region('wen')
  cut = sel.cuts['wen']
  weight = '%s'%(sel.weights['wen']%lumi)
  logger.info('makeLimitForest.py',cut)
  logger.info('makeLimitForest.py',weight)
  processes['wen'] = [
    root.Process('Data',tSingleEle,vmW,dataCut(cut,sel.eleTrigger),'1'),
    root.Process('Wlv',tWlv,vmW,cut,weight),
    root.Process('Zll',tZll,vmW,cut,weight),
    root.Process('ttbar',tTTbar,vmW,cut,weight),
    root.Process('ST',tST,vmW,cut,weight),
    root.Process('Diboson',tVV,vmW,cut,weight),
    root.Process('QCD',tQCD,vmW,cut,weight),
  ]
  processes['wen'] += shiftBtags('ttbar',tTTbar,vmW,cut,'wen')
  processes['wen'] += shiftBtags('Wlv',tWlv,vmW,cut,'wen')
  processes['wen'] += shiftBtags('ST',tST,vmW,cut,'wen')
  processes['wen'] += shiftBtags('Zll',tZll,vmW,cut,'wen')
  processes['wen'] += shiftBtags('Diboson',tVV,vmW,cut,'wen')
  processes['wen'] += shiftBtags('QCD',tQCD,vmW,cut,'wen')
  for p in processes['wen']:
    regions['wen'].AddProcess(p)
  factory.AddRegion(regions['wen'])


# tm
if enable('tm'):
  regions['tm'] = root.Region('tm')
  cut = sel.cuts['tm']
  weight = '%s'%(sel.weights['tm']%lumi)
  logger.info('makeLimitForest.py',cut)
  logger.info('makeLimitForest.py',weight)
  processes['tm'] = [
    root.Process('Data',tMET,vmW,dataCut(cut,sel.metTrigger),'1'),
    root.Process('Wlv',tWlv,vmW,cut,weight),
    root.Process('Zll',tZll,vmW,cut,weight),
    root.Process('ttbar',tTTbar,vmW,cut,weight),
    root.Process('ST',tST,vmW,cut,weight),
    root.Process('Diboson',tVV,vmW,cut,weight),
    root.Process('QCD',tQCD,vmW,cut,weight),
  ]
  processes['tm'] += shiftBtags('ttbar',tTTbar,vmW,cut,'tm')
  processes['tm'] += shiftBtags('Wlv',tWlv,vmW,cut,'tm')
  processes['tm'] += shiftBtags('ST',tST,vmW,cut,'tm')
  processes['tm'] += shiftBtags('Zll',tZll,vmW,cut,'tm')
  processes['tm'] += shiftBtags('Diboson',tVV,vmW,cut,'tm')
  processes['tm'] += shiftBtags('QCD',tQCD,vmW,cut,'tm')
  for p in processes['tm']:
    regions['tm'].AddProcess(p)
  factory.AddRegion(regions['tm'])

# te
if enable('te'):
  regions['te'] = root.Region('te')
  cut = sel.cuts['te']
  weight = '%s'%(sel.weights['te']%lumi)
  logger.info('makeLimitForest.py',cut)
  logger.info('makeLimitForest.py',weight)
  processes['te'] = [
    root.Process('Data',tSingleEle,vmW,dataCut(cut,sel.eleTrigger),'1'),
    root.Process('Wlv',tWlv,vmW,cut,weight),
    root.Process('Zll',tZll,vmW,cut,weight),
    root.Process('ttbar',tTTbar,vmW,cut,weight),
    root.Process('ST',tST,vmW,cut,weight),
    root.Process('Diboson',tVV,vmW,cut,weight),
    root.Process('QCD',tQCD,vmW,cut,weight),
  ]
  processes['te'] += shiftBtags('ttbar',tTTbar,vmW,cut,'te')
  processes['te'] += shiftBtags('Wlv',tWlv,vmW,cut,'te')
  processes['te'] += shiftBtags('ST',tST,vmW,cut,'te')
  processes['te'] += shiftBtags('Zll',tZll,vmW,cut,'te')
  processes['te'] += shiftBtags('Diboson',tVV,vmW,cut,'te')
  processes['te'] += shiftBtags('QCD',tQCD,vmW,cut,'te')
  for p in processes['te']:
    regions['te'].AddProcess(p)
  factory.AddRegion(regions['te'])


# zmm
if enable('zmm'):
  regions['zmm'] = root.Region('zmm')
  cut = sel.cuts['zmm']
  weight = '%s'%(sel.weights['zmm']%lumi)
  logger.info('makeLimitForest.py',cut)
  logger.info('makeLimitForest.py',weight)
  processes['zmm'] = [
    root.Process('Data',tMET,vmZ,dataCut(cut,sel.metTrigger),'1'),
    root.Process('Zvv',tZvv,vmZ,cut,weight),
    root.Process('Wlv',tWlv,vmZ,cut,weight),
    root.Process('Zll',tZll,vmZ,cut,weight),
    root.Process('ttbar',tTTbar,vmZ,cut,weight),
    root.Process('ST',tST,vmZ,cut,weight),
    root.Process('Diboson',tVV,vmZ,cut,weight),
    root.Process('QCD',tQCD,vmZ,cut,weight),
  ]

  for p in processes['zmm']:
    regions['zmm'].AddProcess(p)
  factory.AddRegion(regions['zmm'])

# zee
if enable('zee'):
  regions['zee'] = root.Region('zee')
  cut = sel.cuts['zee']
  weight = '%s'%(sel.weights['zee']%lumi)
  logger.info('makeLimitForest.py',cut)
  logger.info('makeLimitForest.py',weight)
  processes['zee'] = [
    root.Process('Data',tSingleEle,vmZ,dataCut(cut,sel.eleTrigger),'1'),
    root.Process('Zvv',tZvv,vmZ,cut,weight),
    root.Process('Wlv',tWlv,vmZ,cut,weight),
    root.Process('Zll',tZll,vmZ,cut,weight),
    root.Process('ttbar',tTTbar,vmZ,cut,weight),
    root.Process('ST',tST,vmZ,cut,weight),
    root.Process('Diboson',tVV,vmZ,cut,weight),
    root.Process('QCD',tQCD,vmZ,cut,weight),
  ]

  for p in processes['zee']:
    regions['zee'].AddProcess(p)
  factory.AddRegion(regions['zee'])


# photon
if enable('pho'):
  regions['pho'] = root.Region('pho')
  cut = sel.cuts['pho']
  weight = '%s'%(sel.weights['pho']%lumi)
  logger.info('makeLimitForest.py',cut)
  logger.info('makeLimitForest.py',weight)
  processes['pho'] = [
    root.Process('Data',tSinglePho,vmA,dataCut(cut,sel.phoTrigger),'1'),
    root.Process('Pho',tPho,vmA,cut,weight),
    # root.Process('QCD',tSinglePho,vmA,dataCut(cut,phoTrigger),'photonPurityWeight'),
    root.Process('QCD',tQCD,vmA,cut,weight),
  ]
  for p in processes['pho']:
    regions['pho'].AddProcess(p)
  factory.AddRegion(regions['pho'])


logger.info('makeLimitForest','Starting '+str(toProcess))
factory.Run()
logger.info('makeLimitForest','Finishing '+str(toProcess))

factory.Output()

