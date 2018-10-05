from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
}

metFilter='metFilter==1'

#presel = 'jotEta[0]*jotEta[1]<0 && jotPt[0]>80 && jotPt[1]>40 && fabs(jotEta[0])<4.7 && fabs(jotEta[1])<4.7 && (fabs(jotEta[0])<3||fabs(jotEta[0])>3.2) && nTau==0  && jot12Mass>1500 && fabs(jot12DEta)>4.2 && fabs(jot12DPhi)<1.3'
#presel = 'jotEta[0]*jotEta[1]<0 && jotPt[0]>80 && jotPt[1]>40 && fabs(jotEta[0])<4.7 && fabs(jotEta[1])<4.7 && (fabs(jotEta[0])<3||fabs(jotEta[0])>3.2) && nTau==0 && jetNMBtags==0'
#presel = 'jotEta[0]*jotEta[1]<0 && jotPt[0]>80 && jotPt[1]>40 && fabs(jotEta[0])<4.7 && fabs(jotEta[1])<4.7 && nTau==0 && jetNMBtags==0 && jot1VBFID==1'
presel = 'jotEta[0]*jotEta[1]<0 && jotPt[0]>80 && jotPt[1]>40 && fabs(jotEta[0])<4.7 && fabs(jotEta[1])<4.7'
presel = tAND(presel, 'nTau==0 && jetNMBtags==0 && jotVBFID[0]==1 && (fabs(jotEta[0])<3 || fabs(jotEta[1])<3)')
presel = tAND(presel, 'fabs(jotEta[0]+2.82)>0.4 || fabs(jotPhi[0]-2.01)>0.4')
presel = tAND(presel, 'fabs(jotEta[1]+2.82)>0.4 || fabs(jotPhi[1]-2.01)>0.4')
cnc = 'fabs(jot12DEta)>4 && jot12Mass>1300 && fabs(jot12DPhi)<1.5'
mjj = 'fabs(jot12DEta)>1 && fabs(jot12DPhi)<1.5'
nminus1 = 'fabs(jot12DEta)>1'

cuts = {
    'signal'             : tAND(metFilter,tAND(presel,'pfmet>250 && dphipfmet>0.5 && nLooseMuon==0 && nLooseElectron==0 && nLoosePhoton==0 && fabs(calomet-pfmet)/pfmet<0.5')), 
    'singlemuon'         : tAND(metFilter,tAND(presel,'pfUWmag>250 && dphipfUW>0.5 && nLoosePhoton==0 && nLooseMuon==1 && nTightMuon==1 && nLooseElectron==0 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160')),
    'singleelectron'     : tAND(metFilter,tAND(presel,'pfUWmag>250 && dphipfUW>0.5 && nLoosePhoton==0 && nLooseElectron==1 && nTightElectron==1 && nLooseMuon==0 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160 && pfmet>60')),
    'dimuon'             : tAND(metFilter,tAND(presel,'pfUZmag>250 && dphipfUZ>0.5 && nLooseElectron==0 && nLoosePhoton==0 && nLooseMuon==2 && nTightMuon>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5')),
    'dielectron'         : tAND(metFilter,tAND(presel,'pfUZmag>250 && dphipfUZ>0.5 && nLooseMuon==0 && nLoosePhoton==0 && nLooseElectron==2 && nTightElectron>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5')),
    'qcd'                : tAND(metFilter,tAND(presel,'nLooseMuon==0 && nLoosePhoton==0 && nLooseElectron==0')),
}


weights = {
  'signal'         : '%f*sf_pu*sf_l1*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV_VBF*sf_metTrigVBF',
  'w'              : '%f*sf_pu*sf_l1*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV_VBF',
  'z'              : '%f*sf_pu*sf_l1*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV_VBF2l',
  'qcd'            : '%f*sf_pu*sf_l1*normalizedWeight',
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
          weights[x] = tTIMES(weights['z'],'sf_metTrigZmmVBF')
    #      weights[x] = weights['z']
      else:
          weights[x] = tTIMES(weights['w'],'sf_metTrigVBF')

weights_cnc = {k:(v.replace('sf_qcdV_VBF','sf_qcdV_VBFTight')).replace('VBFTight2l','VBF2lTight') 
                for k,v in weights.iteritems()}
