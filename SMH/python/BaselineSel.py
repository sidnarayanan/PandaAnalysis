from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
    'mu':'(trigger&8)!=0',
}

metFilter='metFilter==1'

monojet = 'nJet>0 && jetPt[0]>100'
dijet = 'nJet>1 && jetPt[0]>60 && jetPt[1]>35'

cuts = {
    'singlemu'         : tAND(tAND(metFilter,'dphipfUW>0.5 && nLoosePhoton==0 && nLooseMuon==1 && nTightMuon==1 && nLooseElectron==0 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160'), triggers['mu']),
    'singleel'         : tAND(tAND(metFilter,'dphipfUW>0.5 && nLoosePhoton==0 && nLooseElectron==1 && nTightElectron==1 && nLooseMuon==0 && fabs(calomet-pfmet)/pfUWmag<0.5 && mT<160'), triggers['ele']),
}


