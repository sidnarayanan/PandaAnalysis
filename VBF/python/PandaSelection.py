from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
}

metFilter='metFilter==1 && egmFilter==1'

#presel = 'jot1Eta*jot2Eta<0 && jot1Pt>80 && jot2Pt>40 && fabs(jot1Eta)<4.7 && fabs(jot2Eta)<4.7 && (fabs(jot1Eta)<3||fabs(jot1Eta)>3.2) && nTau==0  && jot12Mass>1500 && fabs(jot12DEta)>4.2 && fabs(jot12DPhi)<1.3'
#presel = 'jot1Eta*jot2Eta<0 && jot1Pt>80 && jot2Pt>40 && fabs(jot1Eta)<4.7 && fabs(jot2Eta)<4.7 && (fabs(jot1Eta)<3||fabs(jot1Eta)>3.2) && nTau==0 && jetNMBtags==0'
#presel = 'jot1Eta*jot2Eta<0 && jot1Pt>80 && jot2Pt>40 && fabs(jot1Eta)<4.7 && fabs(jot2Eta)<4.7 && nTau==0 && jetNMBtags==0 && jot1VBFID==1'
presel = 'jot1Eta*jot2Eta<0 && jot1Pt>80 && jot2Pt>40 && fabs(jot1Eta)<4.7 && fabs(jot2Eta)<4.7 && nTau==0 && jetNMBtags==0 && jot1VBFID==1 && (fabs(jot1Eta)<3 || fabs(jot2Eta)<3)'
cnc = 'fabs(jot12DEta)>4 && jot12Mass>1300 && fabs(jot12DPhi)<1.5'
mjj = 'fabs(jot12DEta)>1 && fabs(jot12DPhi)<1.5'

cuts = {
    'signal'             : tAND(metFilter,tAND(presel,'pfmet>250   && dphipfmet>0.5 && nLooseMuon==0 && nLooseElectron==0 && nLoosePhoton==0 && fabs(calomet-pfmet)/pfmet<0.5')), 
    'singlemuon'         : tAND(metFilter,tAND(presel,'pfUWmag>250 && dphipfUW>0.5 && nLoosePhoton==0 && nLooseMuon==1 && nTightMuon==1 && nLooseElectron==0 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160')),
    'singleelectron'     : tAND(metFilter,tAND(presel,'pfUWmag>250 && dphipfUW>0.5 && nLoosePhoton==0 && nLooseElectron==1 && nTightElectron==1 && nLooseMuon==0 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160 && pfmet>60')),
    'dimuon'             : tAND(metFilter,tAND(presel,'pfUZmag>250 && dphipfUZ>0.5 && nLooseElectron==0 && nLoosePhoton==0 && nLooseMuon==2 && nTightMuon>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5')),
    'dielectron'         : tAND(metFilter,tAND(presel,'pfUZmag>250 && dphipfUZ>0.5 && nLooseMuon==0 && nLoosePhoton==0 && nLooseElectron==2 && nTightElectron>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5')),
    'qcd'                : tAND(metFilter,tAND(presel,'nLooseMuon==0 && nLoosePhoton==0 && nLooseElectron==0')),
}


weights = {
  'signal'         : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV_VBF*sf_metTrigVBF',
  'w'              : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV_VBF',
  'z'              : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV_VBF2l',
  'qcd'            : '%f*sf_pu*normalizedWeight',
}


weights['qcd'] = weights['signal']

for x in ['dimuon','dielectron','singlemuon','singleelectron']:
    if 'electron' in x:
      if 'di' in x:
        weights[x] = weights['z']
      else:
        weights[x] = tTIMES(weights['w'],'sf_eleTrig')
    else:
      if 'di' in x:
          weights[x] = tTIMES(weights['z'],'sf_metTrigZmmVBF')
    #      weights[x] = weights['z']
      else:
          weights[x] = tTIMES(weights['w'],'sf_metTrigVBF')

weights_cnc = {k:(v.replace('sf_qcdV_VBF','sf_qcdV_VBFTight')).replace('VBFTight2l','VBF2lTight') 
                for k,v in weights.iteritems()}
