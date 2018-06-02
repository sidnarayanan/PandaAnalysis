from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
}

metFilter='metFilter==1'
presel = tAND(metFilter,'nFatJet==1 && fj1Pt>250 && fabs(fj1Eta)<2.4 && fj1MSD>50 && nTau==0')

cuts = {
    'tag'           : tAND(presel,'nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && pfUWmag>250 && fj1MaxCSV>0.54 && isojetNBtags==1 && dphipfUW>0.5 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160'),
    'mistag'        : tAND(presel,'nLooseLep==1 && nTightMuon==1 && nLooseElectron==0 && nLoosePhoton==0 && pfUWmag>250 && fj1MaxCSV<0.54 && isojetNBtags==0 && dphipfUW>0.5 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160'),
    'dimuon'        : tAND(presel,'pfUZmag>250 && dphipfUZ>0.5 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && nLooseMuon==2 && nTightLep>0 && 60<diLepMass && diLepMass<120 && fabs(calomet-pfmet)/pfUZmag<0.5'),
    'photon'        : tAND(presel,'nLooseLep==0 && nLoosePhoton==1 && loosePho1IsTight==1 && fabs(loosePho1Eta)<1.4442 && pfUAmag>250'),
}


weights = {
  'tag'            : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_sjbtag1*sf_btag1',
  'mistag'         : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_sjbtag0*sf_btag0',
  'dimuon'         : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV',
  'photon'         : '%f*sf_pu*normalizedWeight*sf_ewkV*sf_qcdV*sf_pho*sf_phoTrig *sf_qcdV2j', # add the additional 2-jet kfactor
}


