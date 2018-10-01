from PandaCore.Tools.Misc import *
from re import sub


triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
}

metFilter='metFilter==1'
topTagSF = '%f*(fj1IsMatched==1)+%f*(fj1IsMatched==0)'%(1.007,1.02)
ak4bTagSF = 'sf_btag0*(isojetNBtags==0)+sf_btag1*(isojetNBtags==1)+1*(isojetNBtags>1)'
photonSF = '0.93'


presel = 'nFatJet==1 && fj1Pt>250 && TMath::Abs(fj1Eta)<2.4 && fj1Tau32<0.61 && 110<fj1MSD && fj1MSD<210'
cuts = {
    # analysis regions
    'signal'            : tAND(metFilter,tAND(presel,'pfmet>175 && puppimet>250 && dphipuppimet>1.1 && (nLooseMuon+nLooseElectron+nLoosePhoton+nTau)==0 && fj1MaxCSV>0.46 && isojetNBtags==0')),
    # 'signal_nomf'       : tAND(presel,'met>175 && puppimet>250 && dphipuppimet>1.1 && (nLooseMuon+nLooseElectron+nLoosePhoton+nTau)==0 && fj1MaxCSV>0.46 && isojetNBtags==0 && fj1isTight==1 && TMath::Abs(met-calomet)/puppimet<0.5'),
    'singlemuontop'     : tAND(metFilter,tAND(presel,'UWmag>250 && (nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==1 && looseLep1IsTight==1 && fj1MaxCSV>0.46 && isojetNBtags==1')),
    'singleelectrontop' : tAND(metFilter,tAND(presel,'UWmag>250 && (nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==1 && looseLep1IsTight==1 && fj1MaxCSV>0.46 && isojetNBtags==1 && puppimet>40')),
    'singlemuonw'       : tAND(metFilter,tAND(presel,'UWmag>250 && (nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==1 && looseLep1IsTight==1 && fj1MaxCSV<0.46 && isojetNBtags==0')),
    'singleelectronw'   : tAND(metFilter,tAND(presel,'UWmag>250 && (nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==1 && looseLep1IsTight==1 && fj1MaxCSV<0.46 && isojetNBtags==0 && puppimet>40')),
    'dimuon'            : tAND(metFilter,tAND(presel,'UZmag>250 && (nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==2 && looseLep1IsTight==1')),
    'dielectron'        : tAND(metFilter,tAND(presel,'UZmag>250 && (nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==2 && looseLep1IsTight==1')),
    'photon'            : tAND(metFilter,tAND(presel,'UAmag>250 && (nLooseMuon+nLooseElectron+nTau)==0 && nLoosePhoton==1 && loosePho1IsTight==1')),
}

tag_presel = removeCut(removeCut(tOR(cuts['singlemuontop'],cuts['singleelectrontop']),'fj1Tau32'),'fj1MSD')
mistag_presel = tAND(removeCut(removeCut(cuts['photon'],'fj1Tau32'),'fj1MSD'),'fj1MSD>40')
tag = 'fj1Tau32<0.61 && 110<fj1MSD && fj1MSD<210'
tt_cuts = {
  'tag' : tag_presel,
  'tag_pass' : tAND(tag,tag_presel),
  'tag_fail' : tAND(tNOT(tag),tag_presel),
  'mistag' : mistag_presel,
  'mistag_pass' : tAND(tag,mistag_presel),
  'mistag_fail' : tAND(tNOT(tag),mistag_presel),
}

'''
  'signal'    : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_ewkZ*sf_qcdZ*sf_ewkW*sf_qcdW*sf_ewkA*sf_qcdA*sf_tt*sf_sjbtag1*sf_lepTrack',topTagSF),ak4bTagSF),
  'top'       : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_ewkZ*sf_qcdZ*sf_ewkW*sf_qcdW*sf_ewkA*sf_qcdA*sf_tt*sf_sjbtag1*sf_lepTrack',topTagSF),ak4bTagSF),
  'w'         : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_ewkZ*sf_qcdZ*sf_ewkW*sf_qcdW*sf_ewkA*sf_qcdA*sf_tt*sf_sjbtag0*sf_lepTrack',topTagSF),ak4bTagSF),
  'notag'     : tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_ewkZ*sf_qcdZ*sf_ewkW*sf_qcdW*sf_ewkA*sf_qcdA*sf_tt*sf_lepTrack',topTagSF),
  'signal_sf' : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_ewkZ*sf_qcdZ*sf_ewkW*sf_qcdW*sf_ewkA*sf_qcdA*sf_tt*sf_sjbtag1',topTagSF),ak4bTagSF),
  'top_sf'    : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_ewkZ*sf_qcdZ*sf_ewkW*sf_qcdW*sf_ewkA*sf_qcdA*sf_tt*sf_sjbtag1',topTagSF),ak4bTagSF),
  'w_sf'      : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_ewkZ*sf_qcdZ*sf_ewkW*sf_qcdW*sf_ewkA*sf_qcdA*sf_tt*sf_sjbtag0',topTagSF),ak4bTagSF),
'''

weights = {
    # analysis weights
  'signal'    : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_ewkV*sf_qcdV*sf_tt*sf_sjbtag1',topTagSF),ak4bTagSF),
  'top'       : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_lepTrack*sf_ewkV*sf_qcdV*sf_tt*sf_sjbtag1',topTagSF),ak4bTagSF),
  'w'         : tTIMES(tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_lepTrack*sf_ewkV*sf_qcdV*sf_tt*sf_sjbtag0',topTagSF),ak4bTagSF),
  'notag'     : tTIMES('%f*normalizedWeight*sf_pu*sf_lep*sf_lepTrack*sf_ewkV*sf_qcdV*sf_tt',topTagSF),
}

for x in ['singlemuontop','singleelectrontop']:
  weights[x] = weights['top']
for x in ['singlemuonw','singleelectronw']:
  weights[x] = weights['w']
for x in ['dimuon','dielectron']:
  weights[x] = weights['notag']
for x in ['photon']:
  weights[x] = tTIMES(photonSF,weights['notag'])

for r in ['signal','top','w','singlemuontop','singleelectrontop','singlemuonw','singleelectronw']:
  for shift in ['BUp','BDown','MUp','MDown']:
    for cent in ['sf_btag','sf_sjbtag']:
      weights[r+'_'+cent+shift] = sub(cent+'0',cent+'0'+shift,sub(cent+'1',cent+'1'+shift,weights[r]))
