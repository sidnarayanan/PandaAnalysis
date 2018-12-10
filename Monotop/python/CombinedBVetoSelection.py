from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
}

metFilter='metFilter==1'
topTagSF = '1'
presel = 'nFatJet==1 && fjPt[0]>250 && fabs(fjEta[0])<2.4'
#presel = 'nFatJet==1 && fjPt[0]>250 && fabs(fjEta[0])<2.4 && 110<fjMSD[0] && fjMSD[0]<210 && 0.1<top_ecf_bdt'

cuts = {
    'signal'          : 'pfmet>250 && dphipfmet>0.5 && nLooseLep==0 && nLoosePhoton==0 && nTau==0 && fabs(calomet-pfmet)/pfmet<0.5 && fjMaxCSV[0]>0.54 && isojetNBtags==0', 
    'singlemuon'      : 'pfUWmag>250 && dphipfUW>0.5 && nLoosePhoton==0 && nTau==0 && nLooseLep==1 && nTightMuon==1 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160',
    'singleelectron'  : 'pfUWmag>250 && dphipfUW>0.5 && nLoosePhoton==0 && nTau==0 && nLooseLep==1 && nTightElectron==1 && pfmet>50 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160',
    'dimuon'          : 'pfUZmag>250 && dphipfUZ>0.5 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && nLooseMuon==2 && nTightMuon>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5 && isojetNBtags==0',
    'dielectron'      : 'pfUZmag>250 && dphipfUZ>0.5 && nLooseMuon==0 && nLoosePhoton==0 && nTau==0 && nLooseElectron==2 && nTightElectron>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5 && isojetNBtags==0',
    'photon'          : 'pfUAmag>250 && dphipfUA>0.5 && nLooseLep==0 && nTau==0 && nLoosePhoton==1 && loosePho1IsTight==1 && fabs(loosePho1Eta)<1.4442 && fabs(calomet-pfmet)/pfUAmag<0.5 && isojetNBtags==0',
}

for k,v in cuts.iteritems():
    cuts[k] = tAND(presel,tAND(metFilter,v))

for r in ['singlemuon','singleelectron']:
	cuts[r+'w'] = tAND(cuts[r],'fjMaxCSV[0]<0.54 && isojetNBtags==0')
	cuts[r+'top'] = tAND(cuts[r],'fjMaxCSV[0]>0.54 && isojetNBtags==1')


weights = {
  'signal'         : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_metTrig*sf_sjbtag1*sf_btag0',
  'top'            : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_sjbtag1*sf_btag1',
  'w'              : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_sjbtag0*sf_btag0',
  'notag'          : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV',
  'z'              : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_btag0',
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
