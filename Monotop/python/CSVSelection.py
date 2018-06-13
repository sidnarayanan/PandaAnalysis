from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
}

metFilter='metFilter==1'

cuts = {
    #'photon'            : tAND(metFilter,'nLoosePhoton==1 && loosePho1IsTight==1 && pfUAmag>200 && (nLooseMuon+nLooseElectron+nTau)==0 && ((nFatJet==1 && fj1MSD<110 && 50<fj1MSD) || nFatJet==0)'),
    'photon'            : tAND(metFilter,'nLoosePhoton==1 && loosePho1IsTight==1 && pfUAmag>200 && (nLooseMuon+nLooseElectron+nTau)==0 && jet1Pt>100 && TMath::Abs(jet1Eta)<2.4'),
    'ttbar'             : tAND(metFilter,'nLooseElectron==1 && nLooseMuon==1 && looseLep1IsTight==1 && (nLoosePhoton+nTau)==0 && pfUWmag>200 && looseLep1PdgId*looseLep2PdgId<0'),
}


weights = {
    'photon'        : '%f*normalizedWeight*sf_pu*sf_ewkV*sf_qcdV*tagCSVWeightv2',
    'ttbar'         : '%f*normalizedWeight*sf_pu*sf_ewkV*sf_qcdV*sf_lep*sf_tt8TeV*fakeCSVWeight*tagCSVWeightv2',
}
