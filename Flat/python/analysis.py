#!/usr/bin/env python

from PandaCore.Utils.load import Load 
from PandaCore.Utils.logging import logger
import ROOT as root

Load('PandaAnalyzer')

def _dump(a):
    logger.info('PandaAnalysis.Flat.analysis','Summary of analysis %s:'%(a.name))
    for k in dir(a):
        if k[0] == '_':
            continue
        if type(getattr(a, k)) != int:
            continue
        logger.info('PandaAnalysis.Flat.analysis','    %20s = %s'%(k, 'True' if getattr(a, k) else 'False'))



def _analysis(name, verbose, **kwargs):
    a = root.pa.Analysis(name)
    for k,v in kwargs.iteritems():
        if not hasattr(a, k):
            logger.error('PandaAnalysis.Flat.analysis','Could not set property %s'%k)
            return None 
        setattr(a, k, bool(v))
    setattr(a, 'dump', lambda : _dump(a))
    if verbose:
        a.dump()
    return a

def analysis(name, verbose=True, **kwargs):
    return _analysis(name, verbose=verbose, **kwargs)


# predefined!
monotop = lambda v=False : _analysis(
        name = 'monotop',
        verbose = v,
    )

dimuon = lambda v=False : _analysis(
        name = 'dimuon',
        verbose = v,
        vbf = True,
        fatjet = False,
        btagSFs = True,
        puppiJets = True,
    )

vbf = lambda v=False : _analysis(
        name = 'vbf',
        verbose = v,
        vbf = True,
        fatjet = False,
        btagSFs = False,
        puppiJets = False
    )

vqqhbb = lambda v=False : _analysis(
        name = 'vbfhbb',
        verbose = v,
        ak8 = True,
        hbb = True,
        vqqhbb = True, 
        fatjet = True,
        btagSFs = True,
        hfCounting = True,
        varyJES = False,
        varyJESTotal = True,
        rerunJES = False,
        rerunJER = True,
        jetFlavorPartons = False,
        jetFlavorJets = True,
    )

vbfhbb = lambda v=False : _analysis(
        name = 'vbfhbb',
        verbose = v,
        ak8 = True,
        hbb = True,
        fatjet = True,
        btagSFs = True,
        btagWeights = True,
        useCMVA = True,
        hfCounting = True,
        bjetBDTReg = False,
        bjetDeepReg = True,
        varyJES = False,
        varyJESTotal = True,
        rerunJES = False,
        rerunJER = True,
        jetFlavorPartons = False,
        jetFlavorJets = True,
        vbf = True,
    )

monoh = lambda v=False : _analysis(
        name = 'monoh',
        verbose = v,
        monoh = True,
    )

gghbb = lambda v=False : _analysis(
        name = 'gghbb',
        verbose = v,
        monoh = True,
        recoil = False,
        ak8 = True,
    )

deep = lambda v=False : _analysis(
        name = 'deep',
        verbose = v,
        ak8 = True,
        deep = True,
        deepTracks = True,
        deepSVs = True,
        deepAntiKtSort = True,
        btagSFs = False,
        jetFlavorPartons = False,
    )

deepgen = lambda v=False : _analysis(
        name = 'deepgen',
        verbose = v,
        ak8 = True,
        btagSFs = False,
        deepGen = True,
        genOnly = True,
    )

kfac = lambda v=False : _analysis(
        name = 'kfac',
        verbose = v,
        hbb = True,
        genOnly = True,
    )

breg = lambda v=False : _analysis(
        name = 'breg',
        verbose = v,
        ak8 = False,
        hbb = True,
        fatjet = False,
        btagSFs = False,
        btagWeights = False,
        useCMVA = False,
        useDeepCSV = True,
        complicatedLeptons = True,
        hfCounting = True,
        bjetRegTraining = True,
        bjetBDTReg = False,
        bjetDeepReg = False,
        varyJES = True,
        rerunJES = True,
        rerunJER = True,
        jetFlavorPartons = False,
        jetFlavorJets = True,
    )

wlnhbb = lambda v=False : _analysis(
        name = 'wlnhbb',
        verbose = v,
        ak8 = True,
        hbb = True,
        monoh = False,
        recoil = True,
        fatjet = True,
        btagSFs = True,
        btagWeights = True,
        useCMVA = False,
        useDeepCSV = True,
        complicatedLeptons = True,
        hfCounting = True,
        recluster = False,
        bjetRegTraining = False,
        bjetBDTReg = False,
        bjetDeepReg = True,
        varyJES = True,
        rerunJES = True,
        rerunJER = True,
        jetFlavorPartons = False,
        jetFlavorJets = True,
        mcTriggers = True,
    )

vv = lambda v=False : _analysis(
        name = 'vv',
        verbose = v,
        ak8 = False,
        hbb = False,
        monoh = False,
        recoil = False,
        fatjet = False,
        btagSFs = True,
        btagWeights = False,
        useCMVA = False,
        useDeepCSV = True,
        complicatedLeptons = True,
        complicatedPhotons = True,
        hfCounting = True,
        recluster = False,
        bjetRegTraining = False,
        bjetBDTReg = True,
        bjetDeepReg = True,
        varyJES = True,
        rerunJES = True,
        rerunJER = True,
        jetFlavorPartons = False,
        jetFlavorJets = True,
        mcTriggers = True,
        vbf = True,
    )
