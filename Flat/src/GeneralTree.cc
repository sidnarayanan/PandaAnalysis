#include "../interface/GeneralTree.h"
// STARTCUSTOM INCLUDE

#include <iostream>

//ENDCUSTOM
GeneralTree::GeneralTree() {
// STARTCUSTOM CONSTRUCTOR
    for (auto ibeta : ibetas) {
        for (auto N : Ns) {
            for (auto order : orders) {
                ECFParams p;
                p.ibeta = ibeta;
                p.N = N;
                p.order = order;
                ecfParams.push_back(p);
                fj1ECFNs[p] = -1;
            }
        }
    }

    for (unsigned iShift=0; iShift!=bNShift; ++iShift) {
        for (unsigned iJet=0; iJet!=bNJet; ++iJet) {
            for (unsigned iTags=0; iTags!=bNTags; ++iTags) {
                BTagParams p;
                p.jet = (BTagJet)iJet;
                p.tag = (BTagTags)iTags;
                p.shift = (BTagShift)iShift;
                btagParams.push_back(p);
                sf_btags[p] = 1;
            }
        }
    }

    for (unsigned iShift=0; iShift<nCsvShifts; iShift++) {
        csvShift shift = csvShifts[iShift];
        sf_csvWeights[shift] = 1;
    }
// ENDCUSTOM
  Reset();
}
~GeneralTree::GeneralTree() {
// STARTCUSTOM DESTRUCTOR
// ENDCUSTOM
}
void GeneralTree::SetAuxTree(TTree *t) {
// STARTCUSTOM AUX
    for (auto p : ecfParams) {
        TString ecfn(makeECFString(p));
        t->Branch("fj1"+ecfn,&(fj1ECFNs[p]),"fj1"+ecfn+"/F");
    }
    t->Branch("fj1Tau32SD",&(fj1Tau32SD),"fj1Tau32SD/F");
    t->Branch("fj1HTTFRec",&(fj1HTTFRec),"fj1HTTFRec/F");
// ENDCUSTOM
}
void GeneralTree::Reset() {
// STARTCUSTOM RESET
    for (auto p : ecfParams) {
        fj1ECFNs[p] = -1;
    }

    for (auto p : btagParams) {
        sf_btags[p] = 1;
    }
    for (unsigned iShift=0; iShift<nCsvShifts; iShift++) {
        csvShift shift = csvShifts[iShift];
        sf_csvWeights[shift] = -1;
    }
    for (auto iter=signal_weights.begin(); iter!=signal_weights.end(); ++iter) {
        signal_weights[iter->first] = 1;
    }
// ENDCUSTOM
  runNumber = 0
  lumiNumber = 0
  eventNumber = -99
  isData = 0
  npv = 0
  pu = 0
  mcWeight = -99
  trigger = 0
  metFilter = 0
  egmFilter = 0
  filter_maxRecoil = -99
  filter_whichRecoil = 0
  badECALFilter = 0
  sf_ewkV = 1
  sf_qcdV = 1
  sf_ewkV2j = 1
  sf_qcdV2j = 1
  sf_qcdV_VBF = 1
  sf_qcdV_VBF2l = 1
  sf_qcdV_VBFTight = 1
  sf_qcdV_VBF2lTight = 1
  sf_qcdTT = 1
  sf_lepID = 1
  sf_lepIso = 1
  sf_lepTrack = 1
  sf_pho = 1
  sf_eleTrig = 1
  sf_muTrig = 1
  sf_phoTrig = 1
  sf_metTrig = 1
  sf_metTrigZmm = 1
  sf_metTrigVBF = 1
  sf_metTrigZmmVBF = 1
  sf_pu = 1
  sf_npv = 1
  sf_tt = 1
  sf_phoPurity = 1
  sumETRaw = -99
  pfmetRaw = -99
  pfmetphi = -99
  pfmetnomu = -99
  puppimet = -99
  puppimetphi = -99
  calomet = -99
  calometphi = -99
  pfcalobalance = -99
  sumET = -99
  trkmet = -99
  trkmetphi = -99
  whichRecoil = 0
  puppiUWphi = -99
  puppiUZphi = -99
  puppiUAphi = -99
  puppiUphi = -99
  pfUWphi = -99
  pfUZphi = -99
  pfUAphi = -99
  pfUphi = -99
  trueGenBosonPt = -99
  genBosonPt = -99
  genBosonEta = -99
  genBosonMass = -99
  genBosonPhi = -99
  genMuonPt = -99
  genMuonEta = -99
  genElectronPt = -99
  genElectronEta = -99
  genTauPt = -99
  genTauEta = -99
  genJet1Pt = -99
  genJet2Pt = -99
  genJet1Eta = -99
  genJet2Eta = -99
  genMjj = -99
  barrelJet1Pt = -99
  barrelJet1Eta = -99
  barrelHT = -99
  barrelHTMiss = -99
  barrelJet12Pt = -99
  jetNBtags = 0
  jetNMBtags = 0
  isojetNBtags = 0
  nFatjet = 0
  fjTau32 = -99
  fjTau21 = -99
  fjTau32SD = -99
  fjTau21SD = -99
  fjMSD = -99
  fjRho = -99
  fjRawRho = -99
  fjRho2 = -99
  fjRawRho2 = -99
  fjMSD_corr = -99
  fjPt = -99
  fjPhi = -99
  fjEta = -99
  fjM = -99
  fjMaxCSV = -99
  fjSubMaxCSV = -99
  fjMinCSV = -99
  fjDoubleCSV = -99
  fjgbb = 0
  fjNbs = 0
  fjGenPt = -99
  fjGenSize = -99
  fjIsMatched = 0
  fjGenWPt = -99
  fjGenWSize = -99
  fjIsWMatched = 0
  fjHighestPtGen = 0
  fjHighestPtGenPt = -99
  fjIsTight = 0
  fjIsLoose = 0
  fjRawPt = -99
  fjNHF = 0
  fjHTTMass = -99
  fjHTTFRec = -99
  fjIsClean = 0
  fjNPartons = 0
  fjNBPartons = 0
  fjNCPartons = 0
  fjPartonM = -99
  fjPartonPt = -99
  fjPartonEta = -99
  nHF = 0
  nB = 0
  nBGenJets = 0
  genFatJetPt = -99
  genFatJetNProngs = 0
  nLoosePhoton = 0
  nTightPhoton = 0
  loosePho1IsTight = 0
  loosePho1Pt = -99
  loosePho1Eta = -99
  loosePho1Phi = -99
  loosePho1SelBit = 0
  looseGenPho1PdgId = 0
  nLooseLep = 0
  nLooseElectron = 0
  nLooseMuon = 0
  nTightLep = 0
  nTightElectron = 0
  nTightMuon = 0
  sf_zz = 1
  sf_zzUnc = 1
  sf_wz = 1
  sf_vh = 1
  sf_vhUp = 1
  sf_vhDown = 1
  genLep1Pt = -99
  genLep1Eta = -99
  genLep1Phi = -99
  genLep1PdgId = 0
  genLep2Pt = -99
  genLep2Eta = -99
  genLep2Phi = -99
  genLep2PdgId = 0
  genLep3Pt = -99
  genLep3Eta = -99
  genLep3Phi = -99
  genLep3PdgId = 0
  genLep4Pt = -99
  genLep4Eta = -99
  genLep4Phi = -99
  genLep4PdgId = 0
  looseGenLep1PdgId = 0
  looseGenLep2PdgId = 0
  looseGenLep3PdgId = 0
  looseGenLep4PdgId = 0
  diLepMass = -99
  nTau = 0
  mT = -99
  hbbeta = -99
  hbbphi = -99
  nSoft2 = 0
  nSoft5 = 0
  nSoft10 = 0
  hbbCosThetaJJ = -99
  hbbCosThetaCSJ1 = -99
  topWBosonCosThetaCS = -99
  topWBosonPt = -99
  topWBosonEta = -99
  topWBosonPhi = -99
  sumEtSoft1 = -99
  scaleUp = -99
  scaleDown = -99
  pdfUp = -99
  pdfDown = -99
  isGS = 0
  for (int iA=0; iA!=NLEP; ++iA) {
    electronPt[iA] = -99
    electronEta[iA] = -99
    electronPhi[iA] = -99
    electronSelBit[iA] = 0
    electronPdgId[iA] = 0
    electronSfLoose[iA] = -99
    electronSfMedium[iA] = -99
    electronSfTight[iA] = -99
    electronSfMvaWP90[iA] = -99
    electronSfMvaWP80[iA] = -99
    electronSfUnc[iA] = -99
    electronSfReco[iA] = -99
    electronD0[iA] = -99
    electronDZ[iA] = -99
    electronNMissingHits[iA] = 0
    electronTripleCharge[iA] = 0
    electronCombIso[iA] = -99
    muonPt[iA] = -99
    muonEta[iA] = -99
    muonPhi[iA] = -99
    muonSelBit[iA] = 0
    muonPdgId[iA] = 0
    muonSfLoose[iA] = -99
    muonSfMedium[iA] = -99
    muonSfTight[iA] = -99
    muonSfUnc[iA] = -99
    muonSfReco[iA] = -99
    muonD0[iA] = -99
    muonDZ[iA] = -99
    muonIsSoftMuon[iA] = 0
    muonCombIso[iA] = -99
  }
  for (int iA=0; iA!=NSUBJET; ++iA) {
    fjsjPt[iA] = -99
    fjsjEta[iA] = -99
    fjsjPhi[iA] = -99
    fjsjM[iA] = -99
    fjsjCSV[iA] = -99
    fjsjQGL[iA] = -99
  }
  for (int iA=0; iA!=2; ++iA) {
    jetPt[iA] = -99
    jetPEta[iA] = -99
    jetPhi[iA] = -99
    jetGenPt[iA] = -99
    jetCSV[iA] = -99
    jetFlav[iA] = 0
    jetIsTight[iA] = 0
    jetIsIso[iA] = 0
    hbbjtidx[iA] = 0
    jotRegFac[iA] = -99
  }
  for (int iA=0; iA!=6; ++iA) {
    scale[iA] = 1
  }
  for (int iA=0; iA!=NJET; ++iA) {
    jotEta[iA] = -99
    jotPhi[iA] = -99
    jotCSV[iA] = -99
    jotVBFID[iA] = 0
    jotCMVA[iA] = -99
    jotIso[iA] = 0
    jotQGL[iA] = -99
    jotLep1Pt[iA] = -99
    jotLep1PtRel[iA] = -99
    jotLep1DeltaR[iA] = -99
    jotTrk1Pt[iA] = -99
    jotVtxPt[iA] = -99
    jotVtxMass[iA] = -99
    jotVtx3DVal[iA] = -99
    jotVtx3DErr[iA] = -99
    jotVtxNtrk[iA] = 0
    jotEMF[iA] = -99
    jotHF[iA] = -99
    jotNLep[iA] = 0
    jotGenPt[iA] = -99
    jotFlav[iA] = 0
  }
  for (int iS=0; iS!=3; ++iS) {
    pfmet[iS] = -99
    puppiUWmag[iS] = -99
    puppiUZmag[iS] = -99
    puppiUAmag[iS] = -99
    puppiUperp[iS] = -99
    puppiUpara[iS] = -99
    pfUWmag[iS] = -99
    pfUZmag[iS] = -99
    puppiUmag[iS] = -99
    pfUAmag[iS] = -99
    pfUperp[iS] = -99
    pfUpara[iS] = -99
    pfUmag[iS] = -99
    dphipfmet[iS] = -99
    dphipuppimet[iS] = -99
    dphipuppiUW[iS] = -99
    dphipuppiUZ[iS] = -99
    dphipuppiUA[iS] = -99
    dphipfUW[iS] = -99
    dphipfUZ[iS] = -99
    dphipfUA[iS] = -99
    dphipuppiU[iS] = -99
    dphipfU[iS] = -99
    nJet[iS] = 0
    nJot[iS] = 0
    nIsoJet[iS] = 0
    jot12Mass[iS] = -99
    jot12DEta[iS] = -99
    jot12DPhi[iS] = -99
    hbbpt[iS] = -99
    hbbm[iS] = -99
    hbbm_reg[iS] = -99
    topMassLep1Met[iS] = -99
  }
  for (int iS=0; iS!=3; ++iS) {
    for (int iA=0; iA!=NJET; ++iA) {
      jotPt[iS][iA] = -99
      jotE[iS][iA] = -99
    }
  }
}
void GeneralTree::WriteTree(TTree *t) {
  treePtr = t;
// STARTCUSTOM WRITE
    for (auto iter=signal_weights.begin(); iter!=signal_weights.end(); ++iter) {
        Book("rw_"+iter->first,&(signal_weights[iter->first]),"rw_"+iter->first+"/F");
    }
    for (auto p : ecfParams) {
        TString ecfn(makeECFString(p));
        Book("fj1"+ecfn,&(fj1ECFNs[p]),"fj1"+ecfn+"/F");
    }

    for (auto p : btagParams) {
        TString btagn(makeBTagSFString(p));
        Book(btagn,&(sf_btags[p]),btagn+"/F");
    }
    if (btagWeights) {
        for (unsigned iShift=0; iShift<nCsvShifts; iShift++) {
            csvShift shift = csvShifts[iShift];
            TString csvWeightString = makeCsvWeightString(shift, useCMVA);
            Book(csvWeightString, &(sf_csvWeights[shift]), csvWeightString+"/F");
        }
    }
// ENDCUSTOM
  Book("runNumber",&runNumber,"runNumber/I");
  Book("lumiNumber",&lumiNumber,"lumiNumber/I");
  Book("eventNumber",&eventNumber,"eventNumber/l");
  Book("isData",&isData,"isData/I");
  Book("npv",&npv,"npv/I");
  Book("pu",&pu,"pu/I");
  Book("mcWeight",&mcWeight,"mcWeight/F");
  Book("trigger",&trigger,"trigger/I");
  Book("metFilter",&metFilter,"metFilter/I");
  Book("egmFilter",&egmFilter,"egmFilter/I");
  Book("filter_maxRecoil",&filter_maxRecoil,"filter_maxRecoil/F");
  Book("filter_whichRecoil",&filter_whichRecoil,"filter_whichRecoil/I");
  Book("badECALFilter",&badECALFilter,"badECALFilter/I");
  Book("sf_ewkV",&sf_ewkV,"sf_ewkV/F");
  Book("sf_qcdV",&sf_qcdV,"sf_qcdV/F");
  Book("sf_ewkV2j",&sf_ewkV2j,"sf_ewkV2j/F");
  Book("sf_qcdV2j",&sf_qcdV2j,"sf_qcdV2j/F");
  Book("sf_qcdV_VBF",&sf_qcdV_VBF,"sf_qcdV_VBF/F");
  Book("sf_qcdV_VBF2l",&sf_qcdV_VBF2l,"sf_qcdV_VBF2l/F");
  Book("sf_qcdV_VBFTight",&sf_qcdV_VBFTight,"sf_qcdV_VBFTight/F");
  Book("sf_qcdV_VBF2lTight",&sf_qcdV_VBF2lTight,"sf_qcdV_VBF2lTight/F");
  Book("sf_qcdTT",&sf_qcdTT,"sf_qcdTT/F");
  Book("sf_lepID",&sf_lepID,"sf_lepID/F");
  Book("sf_lepIso",&sf_lepIso,"sf_lepIso/F");
  Book("sf_lepTrack",&sf_lepTrack,"sf_lepTrack/F");
  Book("sf_pho",&sf_pho,"sf_pho/F");
  Book("sf_eleTrig",&sf_eleTrig,"sf_eleTrig/F");
  Book("sf_muTrig",&sf_muTrig,"sf_muTrig/F");
  Book("sf_phoTrig",&sf_phoTrig,"sf_phoTrig/F");
  Book("sf_metTrig",&sf_metTrig,"sf_metTrig/F");
  Book("sf_metTrigZmm",&sf_metTrigZmm,"sf_metTrigZmm/F");
  Book("sf_metTrigVBF",&sf_metTrigVBF,"sf_metTrigVBF/F");
  Book("sf_metTrigZmmVBF",&sf_metTrigZmmVBF,"sf_metTrigZmmVBF/F");
  Book("sf_pu",&sf_pu,"sf_pu/F");
  Book("sf_npv",&sf_npv,"sf_npv/F");
  Book("sf_tt",&sf_tt,"sf_tt/F");
  Book("sf_phoPurity",&sf_phoPurity,"sf_phoPurity/F");
  Book("sumETRaw",&sumETRaw,"sumETRaw/F");
  Book("pfmetRaw",&pfmetRaw,"pfmetRaw/F");
  Book("pfmet",&(pfmet[0]),"pfmet/F");
  Book("pfmet_JESUp",&(pfmet[1]),"pfmet_JESUp/F");
  Book("pfmet_JESDown",&(pfmet[2]),"pfmet_JESDown/F");
  Book("pfmetphi",&pfmetphi,"pfmetphi/F");
  Book("pfmetnomu",&pfmetnomu,"pfmetnomu/F");
  Book("puppimet",&puppimet,"puppimet/F");
  Book("puppimetphi",&puppimetphi,"puppimetphi/F");
  Book("calomet",&calomet,"calomet/F");
  Book("calometphi",&calometphi,"calometphi/F");
  Book("pfcalobalance",&pfcalobalance,"pfcalobalance/F");
  Book("sumET",&sumET,"sumET/F");
  Book("trkmet",&trkmet,"trkmet/F");
  Book("trkmetphi",&trkmetphi,"trkmetphi/F");
  Book("whichRecoil",&whichRecoil,"whichRecoil/I");
  Book("puppiUWmag",&(puppiUWmag[0]),"puppiUWmag/F");
  Book("puppiUWmag_JESUp",&(puppiUWmag[1]),"puppiUWmag_JESUp/F");
  Book("puppiUWmag_JESDown",&(puppiUWmag[2]),"puppiUWmag_JESDown/F");
  Book("puppiUZmag",&(puppiUZmag[0]),"puppiUZmag/F");
  Book("puppiUZmag_JESUp",&(puppiUZmag[1]),"puppiUZmag_JESUp/F");
  Book("puppiUZmag_JESDown",&(puppiUZmag[2]),"puppiUZmag_JESDown/F");
  Book("puppiUAmag",&(puppiUAmag[0]),"puppiUAmag/F");
  Book("puppiUAmag_JESUp",&(puppiUAmag[1]),"puppiUAmag_JESUp/F");
  Book("puppiUAmag_JESDown",&(puppiUAmag[2]),"puppiUAmag_JESDown/F");
  Book("puppiUperp",&(puppiUperp[0]),"puppiUperp/F");
  Book("puppiUperp_JESUp",&(puppiUperp[1]),"puppiUperp_JESUp/F");
  Book("puppiUperp_JESDown",&(puppiUperp[2]),"puppiUperp_JESDown/F");
  Book("puppiUpara",&(puppiUpara[0]),"puppiUpara/F");
  Book("puppiUpara_JESUp",&(puppiUpara[1]),"puppiUpara_JESUp/F");
  Book("puppiUpara_JESDown",&(puppiUpara[2]),"puppiUpara_JESDown/F");
  Book("pfUWmag",&(pfUWmag[0]),"pfUWmag/F");
  Book("pfUWmag_JESUp",&(pfUWmag[1]),"pfUWmag_JESUp/F");
  Book("pfUWmag_JESDown",&(pfUWmag[2]),"pfUWmag_JESDown/F");
  Book("pfUZmag",&(pfUZmag[0]),"pfUZmag/F");
  Book("pfUZmag_JESUp",&(pfUZmag[1]),"pfUZmag_JESUp/F");
  Book("pfUZmag_JESDown",&(pfUZmag[2]),"pfUZmag_JESDown/F");
  Book("puppiUmag",&(puppiUmag[0]),"puppiUmag/F");
  Book("puppiUmag_JESUp",&(puppiUmag[1]),"puppiUmag_JESUp/F");
  Book("puppiUmag_JESDown",&(puppiUmag[2]),"puppiUmag_JESDown/F");
  Book("pfUAmag",&(pfUAmag[0]),"pfUAmag/F");
  Book("pfUAmag_JESUp",&(pfUAmag[1]),"pfUAmag_JESUp/F");
  Book("pfUAmag_JESDown",&(pfUAmag[2]),"pfUAmag_JESDown/F");
  Book("pfUperp",&(pfUperp[0]),"pfUperp/F");
  Book("pfUperp_JESUp",&(pfUperp[1]),"pfUperp_JESUp/F");
  Book("pfUperp_JESDown",&(pfUperp[2]),"pfUperp_JESDown/F");
  Book("pfUpara",&(pfUpara[0]),"pfUpara/F");
  Book("pfUpara_JESUp",&(pfUpara[1]),"pfUpara_JESUp/F");
  Book("pfUpara_JESDown",&(pfUpara[2]),"pfUpara_JESDown/F");
  Book("pfUmag",&(pfUmag[0]),"pfUmag/F");
  Book("pfUmag_JESUp",&(pfUmag[1]),"pfUmag_JESUp/F");
  Book("pfUmag_JESDown",&(pfUmag[2]),"pfUmag_JESDown/F");
  Book("puppiUWphi",&puppiUWphi,"puppiUWphi/F");
  Book("puppiUZphi",&puppiUZphi,"puppiUZphi/F");
  Book("puppiUAphi",&puppiUAphi,"puppiUAphi/F");
  Book("puppiUphi",&puppiUphi,"puppiUphi/F");
  Book("pfUWphi",&pfUWphi,"pfUWphi/F");
  Book("pfUZphi",&pfUZphi,"pfUZphi/F");
  Book("pfUAphi",&pfUAphi,"pfUAphi/F");
  Book("pfUphi",&pfUphi,"pfUphi/F");
  Book("dphipfmet",&(dphipfmet[0]),"dphipfmet/F");
  Book("dphipfmet_JESUp",&(dphipfmet[1]),"dphipfmet_JESUp/F");
  Book("dphipfmet_JESDown",&(dphipfmet[2]),"dphipfmet_JESDown/F");
  Book("dphipuppimet",&(dphipuppimet[0]),"dphipuppimet/F");
  Book("dphipuppimet_JESUp",&(dphipuppimet[1]),"dphipuppimet_JESUp/F");
  Book("dphipuppimet_JESDown",&(dphipuppimet[2]),"dphipuppimet_JESDown/F");
  Book("dphipuppiUW",&(dphipuppiUW[0]),"dphipuppiUW/F");
  Book("dphipuppiUW_JESUp",&(dphipuppiUW[1]),"dphipuppiUW_JESUp/F");
  Book("dphipuppiUW_JESDown",&(dphipuppiUW[2]),"dphipuppiUW_JESDown/F");
  Book("dphipuppiUZ",&(dphipuppiUZ[0]),"dphipuppiUZ/F");
  Book("dphipuppiUZ_JESUp",&(dphipuppiUZ[1]),"dphipuppiUZ_JESUp/F");
  Book("dphipuppiUZ_JESDown",&(dphipuppiUZ[2]),"dphipuppiUZ_JESDown/F");
  Book("dphipuppiUA",&(dphipuppiUA[0]),"dphipuppiUA/F");
  Book("dphipuppiUA_JESUp",&(dphipuppiUA[1]),"dphipuppiUA_JESUp/F");
  Book("dphipuppiUA_JESDown",&(dphipuppiUA[2]),"dphipuppiUA_JESDown/F");
  Book("dphipfUW",&(dphipfUW[0]),"dphipfUW/F");
  Book("dphipfUW_JESUp",&(dphipfUW[1]),"dphipfUW_JESUp/F");
  Book("dphipfUW_JESDown",&(dphipfUW[2]),"dphipfUW_JESDown/F");
  Book("dphipfUZ",&(dphipfUZ[0]),"dphipfUZ/F");
  Book("dphipfUZ_JESUp",&(dphipfUZ[1]),"dphipfUZ_JESUp/F");
  Book("dphipfUZ_JESDown",&(dphipfUZ[2]),"dphipfUZ_JESDown/F");
  Book("dphipfUA",&(dphipfUA[0]),"dphipfUA/F");
  Book("dphipfUA_JESUp",&(dphipfUA[1]),"dphipfUA_JESUp/F");
  Book("dphipfUA_JESDown",&(dphipfUA[2]),"dphipfUA_JESDown/F");
  Book("dphipuppiU",&(dphipuppiU[0]),"dphipuppiU/F");
  Book("dphipuppiU_JESUp",&(dphipuppiU[1]),"dphipuppiU_JESUp/F");
  Book("dphipuppiU_JESDown",&(dphipuppiU[2]),"dphipuppiU_JESDown/F");
  Book("dphipfU",&(dphipfU[0]),"dphipfU/F");
  Book("dphipfU_JESUp",&(dphipfU[1]),"dphipfU_JESUp/F");
  Book("dphipfU_JESDown",&(dphipfU[2]),"dphipfU_JESDown/F");
  Book("trueGenBosonPt",&trueGenBosonPt,"trueGenBosonPt/F");
  Book("genBosonPt",&genBosonPt,"genBosonPt/F");
  Book("genBosonEta",&genBosonEta,"genBosonEta/F");
  Book("genBosonMass",&genBosonMass,"genBosonMass/F");
  Book("genBosonPhi",&genBosonPhi,"genBosonPhi/F");
  Book("genMuonPt",&genMuonPt,"genMuonPt/F");
  Book("genMuonEta",&genMuonEta,"genMuonEta/F");
  Book("genElectronPt",&genElectronPt,"genElectronPt/F");
  Book("genElectronEta",&genElectronEta,"genElectronEta/F");
  Book("genTauPt",&genTauPt,"genTauPt/F");
  Book("genTauEta",&genTauEta,"genTauEta/F");
  Book("genJet1Pt",&genJet1Pt,"genJet1Pt/F");
  Book("genJet2Pt",&genJet2Pt,"genJet2Pt/F");
  Book("genJet1Eta",&genJet1Eta,"genJet1Eta/F");
  Book("genJet2Eta",&genJet2Eta,"genJet2Eta/F");
  Book("genMjj",&genMjj,"genMjj/F");
  Book("nJet",&(nJet[0]),"nJet/I");
  Book("nJet_JESUp",&(nJet[1]),"nJet_JESUp/I");
  Book("nJet_JESDown",&(nJet[2]),"nJet_JESDown/I");
  Book("nJot",&(nJot[0]),"nJot/I");
  Book("nJot_JESUp",&(nJot[1]),"nJot_JESUp/I");
  Book("nJot_JESDown",&(nJot[2]),"nJot_JESDown/I");
  Book("nIsoJet",&(nIsoJet[0]),"nIsoJet/I");
  Book("nIsoJet_JESUp",&(nIsoJet[1]),"nIsoJet_JESUp/I");
  Book("nIsoJet_JESDown",&(nIsoJet[2]),"nIsoJet_JESDown/I");
  Book("jotPt",jotPt[0],"jotPt[NJET]/F");
  Book("jotPt_JESUp",jotPt[1],"jotPt_JESUp[NJET]/F");
  Book("jotPt_JESDown",jotPt[2],"jotPt_JESDown[NJET]/F");
  Book("jotE",jotE[0],"jotE[NJET]/F");
  Book("jotE_JESUp",jotE[1],"jotE_JESUp[NJET]/F");
  Book("jotE_JESDown",jotE[2],"jotE_JESDown[NJET]/F");
  Book("jotEta",jotEta,"jotEta[NJET]/F");
  Book("jotPhi",jotPhi,"jotPhi[NJET]/F");
  Book("jotCSV",jotCSV,"jotCSV[NJET]/F");
  Book("jotVBFID",jotVBFID,"jotVBFID[NJET]/I");
  Book("barrelJet1Pt",&barrelJet1Pt,"barrelJet1Pt/F");
  Book("barrelJet1Eta",&barrelJet1Eta,"barrelJet1Eta/F");
  Book("barrelHT",&barrelHT,"barrelHT/F");
  Book("barrelHTMiss",&barrelHTMiss,"barrelHTMiss/F");
  Book("barrelJet12Pt",&barrelJet12Pt,"barrelJet12Pt/F");
  Book("jot12Mass",&(jot12Mass[0]),"jot12Mass/F");
  Book("jot12Mass_JESUp",&(jot12Mass[1]),"jot12Mass_JESUp/F");
  Book("jot12Mass_JESDown",&(jot12Mass[2]),"jot12Mass_JESDown/F");
  Book("jot12DEta",&(jot12DEta[0]),"jot12DEta/F");
  Book("jot12DEta_JESUp",&(jot12DEta[1]),"jot12DEta_JESUp/F");
  Book("jot12DEta_JESDown",&(jot12DEta[2]),"jot12DEta_JESDown/F");
  Book("jot12DPhi",&(jot12DPhi[0]),"jot12DPhi/F");
  Book("jot12DPhi_JESUp",&(jot12DPhi[1]),"jot12DPhi_JESUp/F");
  Book("jot12DPhi_JESDown",&(jot12DPhi[2]),"jot12DPhi_JESDown/F");
  Book("jetNBtags",&jetNBtags,"jetNBtags/I");
  Book("jetNMBtags",&jetNMBtags,"jetNMBtags/I");
  Book("isojetNBtags",&isojetNBtags,"isojetNBtags/I");
  Book("nHF",&nHF,"nHF/I");
  Book("nB",&nB,"nB/I");
  Book("nBGenJets",&nBGenJets,"nBGenJets/I");
  Book("genFatJetPt",&genFatJetPt,"genFatJetPt/F");
  Book("genFatJetNProngs",&genFatJetNProngs,"genFatJetNProngs/I");
  Book("nLoosePhoton",&nLoosePhoton,"nLoosePhoton/I");
  Book("nTightPhoton",&nTightPhoton,"nTightPhoton/I");
  Book("loosePho1IsTight",&loosePho1IsTight,"loosePho1IsTight/I");
  Book("loosePho1Pt",&loosePho1Pt,"loosePho1Pt/F");
  Book("loosePho1Eta",&loosePho1Eta,"loosePho1Eta/F");
  Book("loosePho1Phi",&loosePho1Phi,"loosePho1Phi/F");
  Book("loosePho1SelBit",&loosePho1SelBit,"loosePho1SelBit/I");
  Book("looseGenPho1PdgId",&looseGenPho1PdgId,"looseGenPho1PdgId/I");
  Book("nLooseLep",&nLooseLep,"nLooseLep/I");
  Book("nLooseElectron",&nLooseElectron,"nLooseElectron/I");
  Book("nLooseMuon",&nLooseMuon,"nLooseMuon/I");
  Book("nTightLep",&nTightLep,"nTightLep/I");
  Book("nTightElectron",&nTightElectron,"nTightElectron/I");
  Book("nTightMuon",&nTightMuon,"nTightMuon/I");
  Book("electronPt",electronPt,"electronPt[NLEP]/F");
  Book("electronEta",electronEta,"electronEta[NLEP]/F");
  Book("electronPhi",electronPhi,"electronPhi[NLEP]/F");
  Book("electronSelBit",electronSelBit,"electronSelBit[NLEP]/I");
  Book("electronPdgId",electronPdgId,"electronPdgId[NLEP]/I");
  Book("muonPt",muonPt,"muonPt[NLEP]/F");
  Book("muonEta",muonEta,"muonEta[NLEP]/F");
  Book("muonPhi",muonPhi,"muonPhi[NLEP]/F");
  Book("muonSelBit",muonSelBit,"muonSelBit[NLEP]/I");
  Book("muonPdgId",muonPdgId,"muonPdgId[NLEP]/I");
  Book("sf_zz",&sf_zz,"sf_zz/F");
  Book("sf_zzUnc",&sf_zzUnc,"sf_zzUnc/F");
  Book("sf_wz",&sf_wz,"sf_wz/F");
  Book("sf_vh",&sf_vh,"sf_vh/F");
  Book("sf_vhUp",&sf_vhUp,"sf_vhUp/F");
  Book("sf_vhDown",&sf_vhDown,"sf_vhDown/F");
  Book("diLepMass",&diLepMass,"diLepMass/F");
  Book("nTau",&nTau,"nTau/I");
  Book("mT",&mT,"mT/F");
  Book("scaleUp",&scaleUp,"scaleUp/F");
  Book("scaleDown",&scaleDown,"scaleDown/F");
  Book("pdfUp",&pdfUp,"pdfUp/F");
  Book("pdfDown",&pdfDown,"pdfDown/F");
  Book("scale",scale,"scale[6]/F");
  Book("isGS",&isGS,"isGS/I");
  if (is_leptonic) {
    Book("electronSfLoose",electronSfLoose,"electronSfLoose[NLEP]/F");
    Book("electronSfMedium",electronSfMedium,"electronSfMedium[NLEP]/F");
    Book("electronSfTight",electronSfTight,"electronSfTight[NLEP]/F");
    Book("electronSfMvaWP90",electronSfMvaWP90,"electronSfMvaWP90[NLEP]/F");
    Book("electronSfMvaWP80",electronSfMvaWP80,"electronSfMvaWP80[NLEP]/F");
    Book("electronSfUnc",electronSfUnc,"electronSfUnc[NLEP]/F");
    Book("electronSfReco",electronSfReco,"electronSfReco[NLEP]/F");
    Book("electronD0",electronD0,"electronD0[NLEP]/F");
    Book("electronDZ",electronDZ,"electronDZ[NLEP]/F");
    Book("electronNMissingHits",electronNMissingHits,"electronNMissingHits[NLEP]/I");
    Book("electronTripleCharge",electronTripleCharge,"electronTripleCharge[NLEP]/I");
    Book("electronCombIso",electronCombIso,"electronCombIso[NLEP]/F");
    Book("muonSfLoose",muonSfLoose,"muonSfLoose[NLEP]/F");
    Book("muonSfMedium",muonSfMedium,"muonSfMedium[NLEP]/F");
    Book("muonSfTight",muonSfTight,"muonSfTight[NLEP]/F");
    Book("muonSfUnc",muonSfUnc,"muonSfUnc[NLEP]/F");
    Book("muonSfReco",muonSfReco,"muonSfReco[NLEP]/F");
    Book("muonD0",muonD0,"muonD0[NLEP]/F");
    Book("muonDZ",muonDZ,"muonDZ[NLEP]/F");
    Book("muonIsSoftMuon",muonIsSoftMuon,"muonIsSoftMuon[NLEP]/I");
    Book("muonCombIso",muonCombIso,"muonCombIso[NLEP]/F");
    Book("genLep1Pt",&genLep1Pt,"genLep1Pt/F");
    Book("genLep1Eta",&genLep1Eta,"genLep1Eta/F");
    Book("genLep1Phi",&genLep1Phi,"genLep1Phi/F");
    Book("genLep1PdgId",&genLep1PdgId,"genLep1PdgId/I");
    Book("genLep2Pt",&genLep2Pt,"genLep2Pt/F");
    Book("genLep2Eta",&genLep2Eta,"genLep2Eta/F");
    Book("genLep2Phi",&genLep2Phi,"genLep2Phi/F");
    Book("genLep2PdgId",&genLep2PdgId,"genLep2PdgId/I");
    Book("genLep3Pt",&genLep3Pt,"genLep3Pt/F");
    Book("genLep3Eta",&genLep3Eta,"genLep3Eta/F");
    Book("genLep3Phi",&genLep3Phi,"genLep3Phi/F");
    Book("genLep3PdgId",&genLep3PdgId,"genLep3PdgId/I");
    Book("genLep4Pt",&genLep4Pt,"genLep4Pt/F");
    Book("genLep4Eta",&genLep4Eta,"genLep4Eta/F");
    Book("genLep4Phi",&genLep4Phi,"genLep4Phi/F");
    Book("genLep4PdgId",&genLep4PdgId,"genLep4PdgId/I");
    Book("looseGenLep1PdgId",&looseGenLep1PdgId,"looseGenLep1PdgId/I");
    Book("looseGenLep2PdgId",&looseGenLep2PdgId,"looseGenLep2PdgId/I");
    Book("looseGenLep3PdgId",&looseGenLep3PdgId,"looseGenLep3PdgId/I");
    Book("looseGenLep4PdgId",&looseGenLep4PdgId,"looseGenLep4PdgId/I");
  }
  if (is_monotop||is_vbf) {
    Book("jetPt",jetPt,"jetPt[2]/F");
    Book("jetPEta",jetPEta,"jetPEta[2]/F");
    Book("jetPhi",jetPhi,"jetPhi[2]/F");
    Book("jetGenPt",jetGenPt,"jetGenPt[2]/F");
    Book("jetCSV",jetCSV,"jetCSV[2]/F");
    Book("jetFlav",jetFlav,"jetFlav[2]/I");
    Book("jetIsTight",jetIsTight,"jetIsTight[2]/I");
    Book("jetIsIso",jetIsIso,"jetIsIso[2]/I");
  }
  if (is_fatjet) {
    Book("nFatjet",&nFatjet,"nFatjet/I");
    Book("fjTau32",&fjTau32,"fjTau32/F");
    Book("fjTau21",&fjTau21,"fjTau21/F");
    Book("fjTau32SD",&fjTau32SD,"fjTau32SD/F");
    Book("fjTau21SD",&fjTau21SD,"fjTau21SD/F");
    Book("fjMSD",&fjMSD,"fjMSD/F");
    Book("fjRho",&fjRho,"fjRho/F");
    Book("fjRawRho",&fjRawRho,"fjRawRho/F");
    Book("fjRho2",&fjRho2,"fjRho2/F");
    Book("fjRawRho2",&fjRawRho2,"fjRawRho2/F");
    Book("fjMSD_corr",&fjMSD_corr,"fjMSD_corr/F");
    Book("fjPt",&fjPt,"fjPt/F");
    Book("fjPhi",&fjPhi,"fjPhi/F");
    Book("fjEta",&fjEta,"fjEta/F");
    Book("fjM",&fjM,"fjM/F");
    Book("fjMaxCSV",&fjMaxCSV,"fjMaxCSV/F");
    Book("fjSubMaxCSV",&fjSubMaxCSV,"fjSubMaxCSV/F");
    Book("fjMinCSV",&fjMinCSV,"fjMinCSV/F");
    Book("fjDoubleCSV",&fjDoubleCSV,"fjDoubleCSV/F");
    Book("fjgbb",&fjgbb,"fjgbb/I");
    Book("fjNbs",&fjNbs,"fjNbs/I");
    Book("fjGenPt",&fjGenPt,"fjGenPt/F");
    Book("fjGenSize",&fjGenSize,"fjGenSize/F");
    Book("fjIsMatched",&fjIsMatched,"fjIsMatched/I");
    Book("fjGenWPt",&fjGenWPt,"fjGenWPt/F");
    Book("fjGenWSize",&fjGenWSize,"fjGenWSize/F");
    Book("fjIsWMatched",&fjIsWMatched,"fjIsWMatched/I");
    Book("fjHighestPtGen",&fjHighestPtGen,"fjHighestPtGen/I");
    Book("fjHighestPtGenPt",&fjHighestPtGenPt,"fjHighestPtGenPt/F");
    Book("fjIsTight",&fjIsTight,"fjIsTight/I");
    Book("fjIsLoose",&fjIsLoose,"fjIsLoose/I");
    Book("fjRawPt",&fjRawPt,"fjRawPt/F");
    Book("fjNHF",&fjNHF,"fjNHF/I");
    Book("fjHTTMass",&fjHTTMass,"fjHTTMass/F");
    Book("fjHTTFRec",&fjHTTFRec,"fjHTTFRec/F");
    Book("fjIsClean",&fjIsClean,"fjIsClean/I");
    Book("fjNPartons",&fjNPartons,"fjNPartons/I");
    Book("fjNBPartons",&fjNBPartons,"fjNBPartons/I");
    Book("fjNCPartons",&fjNCPartons,"fjNCPartons/I");
    Book("fjPartonM",&fjPartonM,"fjPartonM/F");
    Book("fjPartonPt",&fjPartonPt,"fjPartonPt/F");
    Book("fjPartonEta",&fjPartonEta,"fjPartonEta/F");
  }
  if (is_monohiggs||is_hbb) {
    Book("jotCMVA",jotCMVA,"jotCMVA[NJET]/F");
    Book("jotIso",jotIso,"jotIso[NJET]/I");
    Book("jotQGL",jotQGL,"jotQGL[NJET]/F");
    Book("jotLep1Pt",jotLep1Pt,"jotLep1Pt[NJET]/F");
    Book("jotLep1PtRel",jotLep1PtRel,"jotLep1PtRel[NJET]/F");
    Book("jotLep1DeltaR",jotLep1DeltaR,"jotLep1DeltaR[NJET]/F");
    Book("jotTrk1Pt",jotTrk1Pt,"jotTrk1Pt[NJET]/F");
    Book("jotVtxPt",jotVtxPt,"jotVtxPt[NJET]/F");
    Book("jotVtxMass",jotVtxMass,"jotVtxMass[NJET]/F");
    Book("jotVtx3DVal",jotVtx3DVal,"jotVtx3DVal[NJET]/F");
    Book("jotVtx3DErr",jotVtx3DErr,"jotVtx3DErr[NJET]/F");
    Book("jotVtxNtrk",jotVtxNtrk,"jotVtxNtrk[NJET]/I");
    Book("jotEMF",jotEMF,"jotEMF[NJET]/F");
    Book("jotHF",jotHF,"jotHF[NJET]/F");
    Book("jotNLep",jotNLep,"jotNLep[NJET]/I");
    Book("jotGenPt",jotGenPt,"jotGenPt[NJET]/F");
    Book("jotFlav",jotFlav,"jotFlav[NJET]/I");
    Book("hbbjtidx",hbbjtidx,"hbbjtidx[2]/I");
    Book("jotRegFac",jotRegFac,"jotRegFac[2]/F");
    Book("fjsjPt",fjsjPt,"fjsjPt[NSUBJET]/F");
    Book("fjsjEta",fjsjEta,"fjsjEta[NSUBJET]/F");
    Book("fjsjPhi",fjsjPhi,"fjsjPhi[NSUBJET]/F");
    Book("fjsjM",fjsjM,"fjsjM[NSUBJET]/F");
    Book("fjsjCSV",fjsjCSV,"fjsjCSV[NSUBJET]/F");
    Book("fjsjQGL",fjsjQGL,"fjsjQGL[NSUBJET]/F");
    Book("hbbpt",&(hbbpt[0]),"hbbpt/F");
    Book("hbbpt_JESUp",&(hbbpt[1]),"hbbpt_JESUp/F");
    Book("hbbpt_JESDown",&(hbbpt[2]),"hbbpt_JESDown/F");
    Book("hbbeta",&hbbeta,"hbbeta/F");
    Book("hbbphi",&hbbphi,"hbbphi/F");
    Book("hbbm",&(hbbm[0]),"hbbm/F");
    Book("hbbm_JESUp",&(hbbm[1]),"hbbm_JESUp/F");
    Book("hbbm_JESDown",&(hbbm[2]),"hbbm_JESDown/F");
    Book("hbbm_reg",&(hbbm_reg[0]),"hbbm_reg/F");
    Book("hbbm_reg_JESUp",&(hbbm_reg[1]),"hbbm_reg_JESUp/F");
    Book("hbbm_reg_JESDown",&(hbbm_reg[2]),"hbbm_reg_JESDown/F");
    Book("nSoft2",&nSoft2,"nSoft2/I");
    Book("nSoft5",&nSoft5,"nSoft5/I");
    Book("nSoft10",&nSoft10,"nSoft10/I");
    Book("hbbCosThetaJJ",&hbbCosThetaJJ,"hbbCosThetaJJ/F");
    Book("hbbCosThetaCSJ1",&hbbCosThetaCSJ1,"hbbCosThetaCSJ1/F");
    Book("topMassLep1Met",&(topMassLep1Met[0]),"topMassLep1Met/F");
    Book("topMassLep1Met_JESUp",&(topMassLep1Met[1]),"topMassLep1Met_JESUp/F");
    Book("topMassLep1Met_JESDown",&(topMassLep1Met[2]),"topMassLep1Met_JESDown/F");
    Book("topWBosonCosThetaCS",&topWBosonCosThetaCS,"topWBosonCosThetaCS/F");
    Book("topWBosonPt",&topWBosonPt,"topWBosonPt/F");
    Book("topWBosonEta",&topWBosonEta,"topWBosonEta/F");
    Book("topWBosonPhi",&topWBosonPhi,"topWBosonPhi/F");
    Book("sumEtSoft1",&sumEtSoft1,"sumEtSoft1/F");
  }
}