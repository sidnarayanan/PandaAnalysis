from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
}

metFilter='metFilter==1 && egmFilter==1'

presel = 'jot1Pt>80 && jot1VBFID==1 && fabs(jot1Eta)<2.4'

cuts = {
    'signal'             : tAND(metFilter,tAND(presel,'pfmet>250   && dphipfmet>0.5 && nLooseMuon==0 && nLooseElectron==0 && nLoosePhoton==0 && fabs(calomet-pfmet)/pfmet<0.5')), 
    'singlemuon'         : tAND(metFilter,tAND(presel,'pfUWmag>250 && dphipfUW>0.5 && nLoosePhoton==0 && nLooseMuon==1 && nTightMuon==1 && nLooseElectron==0 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160')),
    'singleelectron'     : tAND(metFilter,tAND(presel,'pfUWmag>250 && dphipfUW>0.5 && nLoosePhoton==0 && nLooseElectron==1 && nTightElectron==1 && nLooseMuon==0 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160 && pfmet>60')),
    'dimuon'             : tAND(metFilter,tAND(presel,'pfUZmag>250 && dphipfUZ>0.5 && nLooseElectron==0 && nLoosePhoton==0 && nLooseMuon==2 && nTightMuon>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5')),
    'dielectron'         : tAND(metFilter,tAND(presel,'pfUZmag>250 && dphipfUZ>0.5 && nLooseMuon==0 && nLoosePhoton==0 && nLooseElectron==2 && nTightElectron>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5')),
    'qcd'                : tAND(metFilter,tAND(presel,'nLooseMuon==0 && nLoosePhoton==0 && nLooseElectron==0')),
}


weights = {
  'signal'         : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_metTrig',
  'w'              : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV',
  'z'              : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV',
  'qcd'            : '%f*sf_pu*normalizedWeight',
}


weights['qcd'] = weights['signal']

for x in ['dimuon','dielectron','singlemuon','singleelectron']:
    if 'electron' in x:
      if 'di' in x:
        weights[x] = tTIMES(weights['z'], 'sf_eleTrig')
      else:
        weights[x] = tTIMES(weights['w'],'sf_eleTrig')
    else:
      if 'di' in x:
          weights[x] = tTIMES(weights['z'],'sf_metTrigZmm')
      else:
          weights[x] = tTIMES(weights['w'],'sf_metTrig')

