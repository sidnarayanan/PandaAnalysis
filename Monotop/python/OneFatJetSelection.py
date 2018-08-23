from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
}

metFilter='metFilter==1 && egmFilter==1'
#metFilter='metFilter==1'
#metFilter = '1==1'
topTagSF = '1'
presel = 'nFatjet==1 && fj1Pt>250 && fabs(fj1Eta)<2.4'

cuts = {
    'signal'             : tAND(metFilter,tAND(presel,'pfmet>250 && dphipfmet>0.5 && nLooseLep==0 && nLoosePhoton==0 && nTau==0 && fabs(calomet-pfmet)/pfmet<0.5 && fj1MaxCSV>0.54 && isojetNBtags==0')), 
    'signal_nob'         : tAND(metFilter,tAND(presel,'pfmet>250 && dphipfmet>0.5 && nLooseLep==0 && nLoosePhoton==0 && nTau==0 && fabs(calomet-pfmet)/pfmet<0.5')), 
    'singlemuon'         : tAND(metFilter,tAND(presel,'pfUWmag>250 && dphipfUW>0.5 && nLoosePhoton==0 && nTau==0 && nLooseLep==1 && looseLep1IsTight==1 && abs(looseLep1PdgId)==13 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160')),
    'singleelectron'     : tAND(metFilter,tAND(presel,'pfUWmag>250 && dphipfUW>0.5 && nLoosePhoton==0 && nTau==0 && nLooseLep==1 && looseLep1IsTight==1 && looseLep1IsHLTSafe==1 && abs(looseLep1PdgId)==11 && pfmet>50 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160')),
    'dimuon_hf'          : tAND(metFilter,tAND(presel,'pfUZmag>250 && dphipfUZ>0.5 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && nLooseMuon==2 && nTightLep>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5 && looseLep1Pt>200 && (fj1MaxCSV>0.54 || isojetNBtags>0)')),
    'dimuon'            : tAND(metFilter,tAND(presel,'pfUZmag>250 && dphipfUZ>0.5 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && nLooseMuon==2 && nTightLep>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5 && isojetNBtags==0')),
    'dielectron'        : tAND(metFilter,tAND(presel,'pfUZmag>250 && dphipfUZ>0.5 && nLooseMuon==0 && nLoosePhoton==0 && nTau==0 && nLooseElectron==2 && nTightLep>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5 && isojetNBtags==0')),
    'photon'            : tAND(metFilter,tAND(presel,'pfUAmag>250 && dphipfUA>0.5 && nLooseLep==0 && nTau==0 && nLoosePhoton==1 && loosePho1IsTight==1 && fabs(loosePho1Eta)<1.4442 && fabs(calomet-pfmet)/pfUAmag<0.5 && isojetNBtags==0')),
    #'qcd'               : tAND(metFilter,tAND(removeCut(removeCut(presel,'top_ecf_bdt'),'fj1MSD'),'pfmet>250 && dphipfmet<0.1 && nLooseLep==0 && nLoosePhoton==0 && nTau==0 && fabs(calomet-pfmet)/pfmet<0.5')), 
    'qcd'               : tAND(metFilter,tAND(removeCut(removeCut(presel,'top_ecf_bdt'),'fj1MSD'),'pfmet>250 && dphipfmet<0.1 && nLooseLep==0 && nLoosePhoton==0 && nTau==0')), 
}
for r in ['singlemuon','singleelectron']:
	cuts[r+'w'] = tAND(cuts[r],'fj1MaxCSV<0.54 && isojetNBtags==0')
	cuts[r+'top'] = tAND(cuts[r],'fj1MaxCSV>0.54 && isojetNBtags==1')


weights = {
  'signal'         : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_metTrig*sf_sjbtag1*sf_btag0',
  'signal_nob'     : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_metTrig',
  'top'            : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_sjbtag1*sf_btag1',
  'w'              : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_sjbtag0*sf_btag0',
  'notag'          : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV',
  'z'              : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_btag0',
  'dimuon_hf'      : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_btag1*sf_sjbtag1',
  'photon'         : '%f*sf_pu*normalizedWeight*sf_ewkV*sf_qcdV*sf_pho*sf_phoTrig *sf_qcdV2j*sf_btag0', # add the additional 2-jet kfactor
}
weights['qcd'] = weights['signal']

for x in ['singlemuontop','singleelectrontop']:
	if 'electron' in x:
	  weights[x] = tTIMES(weights['top'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['top'],'sf_metTrig')
for x in ['singlemuonw','singleelectronw']:
	if 'electron' in x:
	  weights[x] = tTIMES(weights['w'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['w'],'sf_metTrig')
for x in ['dimuon','dielectron']:
	if 'electron' in x:
	  weights[x] = tTIMES(weights['z'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['z'],'sf_metTrig')
for x in ['singlemuon','singleelectron']:
	if 'electron' in x:
	  weights[x] = tTIMES(weights['notag'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['notag'],'sf_metTrig')

for r in ['signal','top','w','singlemuontop','singleelectrontop','singlemuonw','singleelectronw','dimuon','dielectron','photon','z']:
  for shift in ['BUp','BDown','MUp','MDown']:
    for cent in ['sf_btag','sf_sjbtag']:
      weights[r+'_'+cent+shift] = sub(cent+'0',cent+'0'+shift,sub(cent+'1',cent+'1'+shift,weights[r]))
