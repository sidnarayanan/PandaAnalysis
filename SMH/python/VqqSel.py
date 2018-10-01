from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
    'mu':'(trigger&8)!=0',
}

metFilter='metFilter==1'

cuts = {
        'signal' : 'nFatJet>1 && fjPt[0]>300 && fjPt[1]>300 && fjMSD[0]>30 && fjMSD[1]>30'
}

def n2(i):
    return 'fjECFN_2_3_10[{0}] / TMath::Power(fjECFN_1_2_10[{0}], 2)'.format(i)
cuts['pass'] = tAND(cuts['signal'], n2(0)+'<%f && '+n2(1)+'<%f')
cuts['fail'] = tAND(cuts['signal'], n2(0)+'>%f || '+n2(1)+'>%f')
#n2 = 'fjTau21SD[%i]'
#cuts['pass'] = tAND(cuts['signal'], n2%(0)+'<%f && '+n2%(1)+'<%f && fjMSD[fjHiggsIdx]>0.8')
#cuts['fail'] = tAND(cuts['signal'], n2%(0)+'>%f && '+n2%(1)+'>%f')

weights = {
        'signal' : '%f*normalizedWeight'
        }

