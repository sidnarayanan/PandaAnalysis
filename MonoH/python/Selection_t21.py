from PandaCore.Tools.Misc import *
from re import sub


metTrigger='(trigger&1)!=0'
eleTrigger='(trigger&2)!=0'
phoTrigger='(trigger&4)!=0'


cuts = {}
weights = {}
triggers = {}

baseline = 'metFilter==1 && nFatJet==1 && fj1Pt>200 && nTau==0 && fj1Tau21SD<0.55'

#cuts for specific regions
cuts['signal'] = tAND(baseline,'nLooseLep==0 && nLooseElectron==0 && nLoosePhoton==0 && puppimet>200 && pfmet>200 && isojetNBtags==0 && fj1MSD>100 && fj1MSD<150 && dphipuppimet>0.4 && Sum$(jetPt>30 && jetIso)<2 && fj1MaxCSV>0.46 && fj1MinCSV>0.46')
cuts['wmn'] = tAND(baseline,'(nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==1 && looseLep1IsTight==1 && fj1MaxCSV<0.46 && UWmag>200 && isojetNBtags==0 && fj1MSD>50 && fj1MSD<250 && dphipuppimet>0.4')
cuts['wen'] = tAND(baseline,'(nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==1 && looseLep1IsTight==1 && fj1MaxCSV<0.46 && UWmag>200 && isojetNBtags==0 && fj1MSD>50 && fj1MSD<250 && dphipuppimet>0.4 && puppimet>60')
cuts['tm'] = tAND(baseline,'(nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==1 && looseLep1IsTight==1 && fj1MaxCSV>0.46 && UWmag>200 && isojetNBtags==1 && fj1MSD>50 && fj1MSD<250 && dphipuppimet>0.4')
cuts['te'] = tAND(baseline,'(nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==1 && looseLep1IsTight==1 && fj1MaxCSV>0.46 && UWmag>200 && isojetNBtags==1 && fj1MSD>50 && fj1MSD<250 && dphipuppimet>0.4 && puppimet>60')
cuts['zmm'] = tAND(baseline,'(nLooseElectron+nLoosePhoton+nTau)==0 && nLooseMuon==2 && looseLep1IsTight==1 && UZmag>200 && fj1MSD>100 && fj1MSD<150')
cuts['zee'] = tAND(baseline,'(nLooseMuon+nLoosePhoton+nTau)==0 && nLooseElectron==2 && looseLep1IsTight==1 && UZmag>200 && fj1MSD>100 && fj1MSD<150')
cuts['pho'] = tAND(baseline,'(nLooseMuon+nLooseElectron+nTau)==0 && nLoosePhoton==1 && loosePho1IsTight==1 && UAmag>200 && fj1MSD>100 && fj1MSD<150')

#weights for specific regions
weights['base'] = 'normalizedWeight*sf_pu*sf_ewkV*sf_qcdV*sf_tt*sf_lep*sf_lepTrack'
weights['signal'] = tTIMES(weights['base'],'%f*sf_sjbtag2*sf_btag0*sf_metTrig')
weights['wmn'] = tTIMES(weights['base'],'%f*sf_sjbtag0*sf_btag0*sf_metTrig')
weights['wen'] = tTIMES(weights['base'],'%f*sf_sjbtag0*sf_btag0*sf_eleTrig')
weights['zmm'] = tTIMES(weights['base'],'%f*sf_metTrig')
weights['zee'] = tTIMES(weights['base'],'%f*sf_eleTrig')
weights['pho'] = tTIMES(weights['base'],'%f*sf_phoTrig*0.93')
weights['tm'] = tTIMES(weights['base'],'%f*sf_sjbtag1*sf_btag1*sf_metTrig')
weights['te'] = tTIMES(weights['base'],'%f*sf_sjbtag1*sf_btag1*sf_eleTrig')

for r in ['signal','wmn','wen','tm','te']:
    print r
    for shift in ['BUp','BDown','MUp','MDown']:
        for cent in ['sf_btag','sf_sjbtag']:
            if '2' in weights[r]:
                weights[r+'_'+cent+shift] = sub(cent+'0',cent+'0'+shift,sub(cent+'2',cent+'2'+shift,weights[r]))
#                print weights[r+'_'+cent+shift]
            if '1' in weights[r] or r=='wmn' or r=='wen':
                weights[r+'_'+cent+shift] = sub(cent+'0',cent+'0'+shift,sub(cent+'1',cent+'1'+shift,weights[r]))
#                print weights[r+'_'+cent+shift]
            
