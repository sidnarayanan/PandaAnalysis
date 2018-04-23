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
                fjECFNs[p] = -1;
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
GeneralTree::~GeneralTree() {
// STARTCUSTOM DESTRUCTOR
// ENDCUSTOM
}
void GeneralTree::SetAuxTree(TTree *t) {
// STARTCUSTOM AUX
    for (auto p : ecfParams) {
        TString ecfn(makeECFString(p));
        t->Branch("fj"+ecfn,&(fjECFNs[p]),"fj"+ecfn+"/F");
    }
    t->Branch("fjTau32SD",&(fjTau32SD),"fjTau32SD/F");
    t->Branch("fjHTTFRec",&(fjHTTFRec),"fjHTTFRec/F");
// ENDCUSTOM
}
void GeneralTree::Reset() {
// STARTCUSTOM RESET
    for (auto p : ecfParams) {
        fjECFNs[p] = -1;
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
  runNumber = 0;
  lumiNumber = 0;
  eventNumber = -99;
  isData = 0;
  npv = 0;
  pu = 0;
  mcWeight = -99;
  trigger = 0;
  metFilter = 0;
  egmFilter = 0;
  filter_maxRecoil = -99;
  filter_whichRecoil = 0;
  badECALFilter = 1;
  sf_ewkV = 1;
  sf_qcdV = 1;
  sf_ewkV2j = 1;
  sf_qcdV2j = 1;
  sf_qcdV_VBF = 1;
  sf_qcdV_VBF2l = 1;
  sf_qcdV_VBFTight = 1;
  sf_qcdV_VBF2lTight = 1;
  sf_qcdTT = 1;
  sf_lepID = 1;
  sf_lepIso = 1;
  sf_lepTrack = 1;
  sf_pho = 1;
  sf_eleTrig = 1;
  sf_muTrig = 1;
  sf_phoTrig = 1;
  sf_metTrig = 1;
  sf_metTrigZmm = 1;
  sf_metTrigVBF = 1;
  sf_metTrigZmmVBF = 1;
  sf_pu = 1;
  sf_npv = 1;
  sf_tt = 1;
  sf_phoPurity = 1;
  sumETRaw = -99;
  pfmetRaw = -99;
  calomet = -99;
  calometphi = -99;
  pfcalobalance = -99;
  sumET = -99;
  trkmet = -99;
  trkmetphi = -99;
  pfmetsig = -99;
  puppimetsig = -99;
  whichRecoil = 0;
  trueGenBosonPt = -99;
  genBosonPt = -99;
  genBosonEta = -99;
  genBosonMass = -99;
  genBosonPhi = -99;
  genMuonPt = -99;
  genMuonEta = -99;
  genElectronPt = -99;
  genElectronEta = -99;
  genTauPt = -99;
  genTauEta = -99;
  genJet1Pt = -99;
  genJet2Pt = -99;
  genJet1Eta = -99;
  genJet2Eta = -99;
  genMjj = -99;
  genTopPt = -99;
  genAntiTopPt = -99;
  nJotMax = 0;
  barrelJet1Pt = -99;
  barrelJet1Eta = -99;
  barrelHT = 0;
  barrelHTMiss = -99;
  barrelJet12Pt = 0;
  nFatjet = 0;
  fjTau32 = -99;
  fjTau21 = -99;
  fjTau32SD = -99;
  fjTau21SD = -99;
  fjMSD = -99;
  fjRho = -99;
  fjRawRho = -99;
  fjRho2 = -99;
  fjRawRho2 = -99;
  fjMSD_corr = -99;
  fjPt = -99;
  fjPhi = -99;
  fjEta = -99;
  fjM = -99;
  fjMaxCSV = -99;
  fjSubMaxCSV = -99;
  fjMinCSV = -99;
  fjDoubleCSV = -99;
  fjgbb = 0;
  fjNbs = 0;
  fjGenPt = -99;
  fjGenSize = -99;
  fjIsMatched = 0;
  fjGenWPt = -99;
  fjGenWSize = -99;
  fjIsWMatched = 0;
  fjHighestPtGen = 0;
  fjHighestPtGenPt = -99;
  fjIsTight = 0;
  fjIsLoose = 0;
  fjRawPt = -99;
  fjNHF = 0;
  fjHTTMass = -99;
  fjHTTFRec = -99;
  fjIsClean = 0;
  fjNPartons = 0;
  fjNBPartons = 0;
  fjNCPartons = 0;
  fjPartonM = -99;
  fjPartonPt = -99;
  fjPartonEta = -99;
  fjGenNumB = 0;
  nHF = 0;
  nB = 0;
  nBGenJets = 0;
  genFatJetPt = -99;
  genFatJetNProngs = 0;
  nLoosePhoton = 0;
  nTightPhoton = 0;
  loosePho1IsTight = 0;
  loosePho1Pt = -99;
  loosePho1Eta = -99;
  loosePho1Phi = -99;
  loosePho1SelBit = 0;
  looseGenPho1PdgId = 0;
  nLooseLep = 0;
  nLooseElectron = 0;
  nLooseMuon = 0;
  nTightLep = 0;
  nTightElectron = 0;
  nTightMuon = 0;
  sf_zz = 1;
  sf_zzUnc = 1;
  sf_wz = 1;
  sf_vh = 1;
  sf_vhUp = 1;
  sf_vhDown = 1;
  genLep1Pt = -99;
  genLep1Eta = -99;
  genLep1Phi = -99;
  genLep1PdgId = 0;
  genLep2Pt = -99;
  genLep2Eta = -99;
  genLep2Phi = -99;
  genLep2PdgId = 0;
  genLep3Pt = -99;
  genLep3Eta = -99;
  genLep3Phi = -99;
  genLep3PdgId = 0;
  genLep4Pt = -99;
  genLep4Eta = -99;
  genLep4Phi = -99;
  genLep4PdgId = 0;
  genWPlusPt = -99;
  genWMinusPt = -99;
  genWPlusEta = -99;
  genWMinusEta = -99;
  looseGenLep1PdgId = 0;
  looseGenLep2PdgId = 0;
  looseGenLep3PdgId = 0;
  looseGenLep4PdgId = 0;
  diLepMass = -99;
  nTau = 0;
  nSoft2 = 0;
  nSoft5 = 0;
  nSoft10 = 0;
  topWBosonPt = -99;
  topWBosonEta = -99;
  topWBosonPhi = -99;
  sumEtSoft1 = -99;
  scaleUp = 1;
  scaleDown = 1;
  pdfUp = 1;
  pdfDown = 1;
  lheHT = -99;
  isGS = 0;
  for (int iA=0; iA!=4; ++iA) {
    electronPt[iA] = -99;
    electronEta[iA] = -99;
    electronPhi[iA] = -99;
    electronSelBit[iA] = 0;
    electronPdgId[iA] = 0;
    electronSfLoose[iA] = -99;
    electronSfMedium[iA] = -99;
    electronSfTight[iA] = -99;
    electronSfMvaWP90[iA] = -99;
    electronSfMvaWP80[iA] = -99;
    electronSfUnc[iA] = -99;
    electronSfReco[iA] = -99;
    electronD0[iA] = -99;
    electronDZ[iA] = -99;
    electronNMissingHits[iA] = 0;
    electronTripleCharge[iA] = 0;
    electronCombIso[iA] = -99;
    muonPt[iA] = -99;
    muonEta[iA] = -99;
    muonPhi[iA] = -99;
    muonSelBit[iA] = 0;
    muonPdgId[iA] = 0;
    muonSfLoose[iA] = -99;
    muonSfMedium[iA] = -99;
    muonSfTight[iA] = -99;
    muonSfUnc[iA] = -99;
    muonSfReco[iA] = -99;
    muonD0[iA] = -99;
    muonDZ[iA] = -99;
    muonIsSoftMuon[iA] = 0;
    muonCombIso[iA] = -99;
  }
  for (int iA=0; iA!=2; ++iA) {
    jetPt[iA] = -99;
    jetEta[iA] = -99;
    jetPhi[iA] = -99;
    jetGenPt[iA] = -99;
    jetCSV[iA] = -99;
    jetFlav[iA] = 0;
    jetIsTight[iA] = 0;
    jetIsIso[iA] = 0;
  }
  for (int iA=0; iA!=2; ++iA) {
    fjsjPt[iA] = -99;
    fjsjEta[iA] = -99;
    fjsjPhi[iA] = -99;
    fjsjM[iA] = -99;
    fjsjCSV[iA] = -99;
    fjsjQGL[iA] = -99;
  }
  for (int iA=0; iA!=6; ++iA) {
    scale[iA] = 1;
  }
  for (int iS=0; iS!=3; ++iS) {
    pfmet[iS] = -99;
    pfmetphi[iS] = -99;
    pfmetnomu[iS] = -99;
    puppimet[iS] = -99;
    puppimetphi[iS] = -99;
    puppiUWmag[iS] = -99;
    puppiUZmag[iS] = -99;
    puppiUAmag[iS] = -99;
    puppiUperp[iS] = -99;
    puppiUpara[iS] = -99;
    puppiUmag[iS] = -99;
    pfUWmag[iS] = -99;
    pfUZmag[iS] = -99;
    pfUAmag[iS] = -99;
    pfUperp[iS] = -99;
    pfUpara[iS] = -99;
    pfUmag[iS] = -99;
    puppiUWphi[iS] = -99;
    puppiUZphi[iS] = -99;
    puppiUAphi[iS] = -99;
    puppiUphi[iS] = -99;
    pfUWphi[iS] = -99;
    pfUZphi[iS] = -99;
    pfUAphi[iS] = -99;
    pfUphi[iS] = -99;
    dphipfmet[iS] = 999;
    dphipuppimet[iS] = 999;
    dphipuppiUW[iS] = 999;
    dphipuppiUZ[iS] = 999;
    dphipuppiUA[iS] = 999;
    dphipfUW[iS] = 999;
    dphipfUZ[iS] = 999;
    dphipfUA[iS] = 999;
    dphipuppiU[iS] = 999;
    dphipfU[iS] = 999;
    nJet[iS] = 0;
    nJot[iS] = 0;
    nIsoJet[iS] = 0;
    jot12Mass[iS] = -99;
    jot12DEta[iS] = -99;
    jot12DPhi[iS] = -99;
    jetNBtags[iS] = 0;
    jetNMBtags[iS] = 0;
    isojetNBtags[iS] = 0;
    mT[iS] = -99;
    hbbpt[iS] = -99;
    hbbeta[iS] = -99;
    hbbphi[iS] = -99;
    hbbm[iS] = -99;
    hbbm_reg[iS] = -99;
    hbbpt_reg[iS] = -99;
    hbbCosThetaJJ[iS] = -99;
    hbbCosThetaCSJ1[iS] = -99;
    topMassLep1Met[iS] = -99;
    topWBosonCosThetaCS[iS] = -99;
  }
  for (int iS=0; iS!=3; ++iS) {
    for (int iA=0; iA!=20; ++iA) {
      jotM[iS][iA] = -99;
      jotCMVA[iS][iA] = -99;
      jotIso[iS][iA] = 0;
      jotQGL[iS][iA] = -99;
      jotLep1Pt[iS][iA] = -99;
      jotLep1PtRel[iS][iA] = -99;
      jotLep1DeltaR[iS][iA] = -99;
      jotTrk1Pt[iS][iA] = -99;
      jotVtxPt[iS][iA] = -99;
      jotVtxMass[iS][iA] = -99;
      jotVtx3DVal[iS][iA] = -99;
      jotVtx3DErr[iS][iA] = -99;
      jotVtxNtrk[iS][iA] = 0;
      jotEMF[iS][iA] = -99;
      jotHF[iS][iA] = -99;
      jotNLep[iS][iA] = 0;
      jotGenPt[iS][iA] = -99;
      jotFlav[iS][iA] = 0;
    }
  }
  for (int iS=0; iS!=3; ++iS) {
    for (int iA=0; iA!=20; ++iA) {
      jotPt[iS][iA] = -99;
      jotE[iS][iA] = -99;
      jotEta[iS][iA] = -99;
      jotPhi[iS][iA] = -99;
      jotCSV[iS][iA] = -99;
      jotVBFID[iS][iA] = 0;
    }
  }
  for (int iS=0; iS!=3; ++iS) {
    for (int iA=0; iA!=2; ++iA) {
      hbbjtidx[iS][iA] = 0;
      hbbjtRegFac[iS][iA] = -99;
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
        Book("fj"+ecfn,&(fjECFNs[p]),"fj"+ecfn+"/F");
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
  Book("pfmetphi",&(pfmetphi[0]),"pfmetphi/F");
  Book("pfmetphi_JESUp",&(pfmetphi[1]),"pfmetphi_JESUp/F");
  Book("pfmetphi_JESDown",&(pfmetphi[2]),"pfmetphi_JESDown/F");
  Book("pfmetnomu",&(pfmetnomu[0]),"pfmetnomu/F");
  Book("pfmetnomu_JESUp",&(pfmetnomu[1]),"pfmetnomu_JESUp/F");
  Book("pfmetnomu_JESDown",&(pfmetnomu[2]),"pfmetnomu_JESDown/F");
  Book("puppimet",&(puppimet[0]),"puppimet/F");
  Book("puppimet_JESUp",&(puppimet[1]),"puppimet_JESUp/F");
  Book("puppimet_JESDown",&(puppimet[2]),"puppimet_JESDown/F");
  Book("puppimetphi",&(puppimetphi[0]),"puppimetphi/F");
  Book("puppimetphi_JESUp",&(puppimetphi[1]),"puppimetphi_JESUp/F");
  Book("puppimetphi_JESDown",&(puppimetphi[2]),"puppimetphi_JESDown/F");
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
  Book("puppiUmag",&(puppiUmag[0]),"puppiUmag/F");
  Book("puppiUmag_JESUp",&(puppiUmag[1]),"puppiUmag_JESUp/F");
  Book("puppiUmag_JESDown",&(puppiUmag[2]),"puppiUmag_JESDown/F");
  Book("pfUWmag",&(pfUWmag[0]),"pfUWmag/F");
  Book("pfUWmag_JESUp",&(pfUWmag[1]),"pfUWmag_JESUp/F");
  Book("pfUWmag_JESDown",&(pfUWmag[2]),"pfUWmag_JESDown/F");
  Book("pfUZmag",&(pfUZmag[0]),"pfUZmag/F");
  Book("pfUZmag_JESUp",&(pfUZmag[1]),"pfUZmag_JESUp/F");
  Book("pfUZmag_JESDown",&(pfUZmag[2]),"pfUZmag_JESDown/F");
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
  Book("puppiUWphi",&(puppiUWphi[0]),"puppiUWphi/F");
  Book("puppiUWphi_JESUp",&(puppiUWphi[1]),"puppiUWphi_JESUp/F");
  Book("puppiUWphi_JESDown",&(puppiUWphi[2]),"puppiUWphi_JESDown/F");
  Book("puppiUZphi",&(puppiUZphi[0]),"puppiUZphi/F");
  Book("puppiUZphi_JESUp",&(puppiUZphi[1]),"puppiUZphi_JESUp/F");
  Book("puppiUZphi_JESDown",&(puppiUZphi[2]),"puppiUZphi_JESDown/F");
  Book("puppiUAphi",&(puppiUAphi[0]),"puppiUAphi/F");
  Book("puppiUAphi_JESUp",&(puppiUAphi[1]),"puppiUAphi_JESUp/F");
  Book("puppiUAphi_JESDown",&(puppiUAphi[2]),"puppiUAphi_JESDown/F");
  Book("puppiUphi",&(puppiUphi[0]),"puppiUphi/F");
  Book("puppiUphi_JESUp",&(puppiUphi[1]),"puppiUphi_JESUp/F");
  Book("puppiUphi_JESDown",&(puppiUphi[2]),"puppiUphi_JESDown/F");
  Book("pfUWphi",&(pfUWphi[0]),"pfUWphi/F");
  Book("pfUWphi_JESUp",&(pfUWphi[1]),"pfUWphi_JESUp/F");
  Book("pfUWphi_JESDown",&(pfUWphi[2]),"pfUWphi_JESDown/F");
  Book("pfUZphi",&(pfUZphi[0]),"pfUZphi/F");
  Book("pfUZphi_JESUp",&(pfUZphi[1]),"pfUZphi_JESUp/F");
  Book("pfUZphi_JESDown",&(pfUZphi[2]),"pfUZphi_JESDown/F");
  Book("pfUAphi",&(pfUAphi[0]),"pfUAphi/F");
  Book("pfUAphi_JESUp",&(pfUAphi[1]),"pfUAphi_JESUp/F");
  Book("pfUAphi_JESDown",&(pfUAphi[2]),"pfUAphi_JESDown/F");
  Book("pfUphi",&(pfUphi[0]),"pfUphi/F");
  Book("pfUphi_JESUp",&(pfUphi[1]),"pfUphi_JESUp/F");
  Book("pfUphi_JESDown",&(pfUphi[2]),"pfUphi_JESDown/F");
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
  Book("genTopPt",&genTopPt,"genTopPt/F");
  Book("genAntiTopPt",&genAntiTopPt,"genAntiTopPt/F");
  Book("nJet",&(nJet[0]),"nJet/I");
  Book("nJet_JESUp",&(nJet[1]),"nJet_JESUp/I");
  Book("nJet_JESDown",&(nJet[2]),"nJet_JESDown/I");
  Book("nJot",&(nJot[0]),"nJot/I");
  Book("nJot_JESUp",&(nJot[1]),"nJot_JESUp/I");
  Book("nJot_JESDown",&(nJot[2]),"nJot_JESDown/I");
  Book("nIsoJet",&(nIsoJet[0]),"nIsoJet/I");
  Book("nIsoJet_JESUp",&(nIsoJet[1]),"nIsoJet_JESUp/I");
  Book("nIsoJet_JESDown",&(nIsoJet[2]),"nIsoJet_JESDown/I");
  Book("jotPt",jotPt[0],"jotPt["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotPt_JESUp",jotPt[1],"jotPt_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotPt_JESDown",jotPt[2],"jotPt_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotE",jotE[0],"jotE["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotE_JESUp",jotE[1],"jotE_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotE_JESDown",jotE[2],"jotE_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotEta",jotEta[0],"jotEta["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotEta_JESUp",jotEta[1],"jotEta_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotEta_JESDown",jotEta[2],"jotEta_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotPhi",jotPhi[0],"jotPhi["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotPhi_JESUp",jotPhi[1],"jotPhi_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotPhi_JESDown",jotPhi[2],"jotPhi_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotCSV",jotCSV[0],"jotCSV["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotCSV_JESUp",jotCSV[1],"jotCSV_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotCSV_JESDown",jotCSV[2],"jotCSV_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
  Book("jotVBFID",jotVBFID[0],"jotVBFID["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
  Book("jotVBFID_JESUp",jotVBFID[1],"jotVBFID_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
  Book("jotVBFID_JESDown",jotVBFID[2],"jotVBFID_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
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
  Book("jetNBtags",&(jetNBtags[0]),"jetNBtags/I");
  Book("jetNBtags_JESUp",&(jetNBtags[1]),"jetNBtags_JESUp/I");
  Book("jetNBtags_JESDown",&(jetNBtags[2]),"jetNBtags_JESDown/I");
  Book("jetNMBtags",&(jetNMBtags[0]),"jetNMBtags/I");
  Book("jetNMBtags_JESUp",&(jetNMBtags[1]),"jetNMBtags_JESUp/I");
  Book("jetNMBtags_JESDown",&(jetNMBtags[2]),"jetNMBtags_JESDown/I");
  Book("isojetNBtags",&(isojetNBtags[0]),"isojetNBtags/I");
  Book("isojetNBtags_JESUp",&(isojetNBtags[1]),"isojetNBtags_JESUp/I");
  Book("isojetNBtags_JESDown",&(isojetNBtags[2]),"isojetNBtags_JESDown/I");
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
  Book("electronPt",electronPt,"electronPt["+TString("nLooseElectron")+"]/F");
  Book("electronEta",electronEta,"electronEta["+TString("nLooseElectron")+"]/F");
  Book("electronPhi",electronPhi,"electronPhi["+TString("nLooseElectron")+"]/F");
  Book("electronSelBit",electronSelBit,"electronSelBit["+TString("nLooseElectron")+"]/I");
  Book("electronPdgId",electronPdgId,"electronPdgId["+TString("nLooseElectron")+"]/I");
  Book("muonPt",muonPt,"muonPt["+TString("nLooseMuon")+"]/F");
  Book("muonEta",muonEta,"muonEta["+TString("nLooseMuon")+"]/F");
  Book("muonPhi",muonPhi,"muonPhi["+TString("nLooseMuon")+"]/F");
  Book("muonSelBit",muonSelBit,"muonSelBit["+TString("nLooseMuon")+"]/I");
  Book("muonPdgId",muonPdgId,"muonPdgId["+TString("nLooseMuon")+"]/I");
  Book("sf_zz",&sf_zz,"sf_zz/F");
  Book("sf_zzUnc",&sf_zzUnc,"sf_zzUnc/F");
  Book("sf_wz",&sf_wz,"sf_wz/F");
  Book("sf_vh",&sf_vh,"sf_vh/F");
  Book("sf_vhUp",&sf_vhUp,"sf_vhUp/F");
  Book("sf_vhDown",&sf_vhDown,"sf_vhDown/F");
  Book("diLepMass",&diLepMass,"diLepMass/F");
  Book("nTau",&nTau,"nTau/I");
  Book("mT",&(mT[0]),"mT/F");
  Book("mT_JESUp",&(mT[1]),"mT_JESUp/F");
  Book("mT_JESDown",&(mT[2]),"mT_JESDown/F");
  Book("scaleUp",&scaleUp,"scaleUp/F");
  Book("scaleDown",&scaleDown,"scaleDown/F");
  Book("pdfUp",&pdfUp,"pdfUp/F");
  Book("pdfDown",&pdfDown,"pdfDown/F");
  Book("scale",scale,"scale["+TString("6")+"]/F");
  Book("lheHT",&lheHT,"lheHT/F");
  Book("isGS",&isGS,"isGS/I");
  if (is_leptonic) {
    Book("electronSfLoose",electronSfLoose,"electronSfLoose["+TString("nLooseElectron")+"]/F");
    Book("electronSfMedium",electronSfMedium,"electronSfMedium["+TString("nLooseElectron")+"]/F");
    Book("electronSfTight",electronSfTight,"electronSfTight["+TString("nLooseElectron")+"]/F");
    Book("electronSfMvaWP90",electronSfMvaWP90,"electronSfMvaWP90["+TString("nLooseElectron")+"]/F");
    Book("electronSfMvaWP80",electronSfMvaWP80,"electronSfMvaWP80["+TString("nLooseElectron")+"]/F");
    Book("electronSfUnc",electronSfUnc,"electronSfUnc["+TString("nLooseElectron")+"]/F");
    Book("electronSfReco",electronSfReco,"electronSfReco["+TString("nLooseElectron")+"]/F");
    Book("electronD0",electronD0,"electronD0["+TString("nLooseElectron")+"]/F");
    Book("electronDZ",electronDZ,"electronDZ["+TString("nLooseElectron")+"]/F");
    Book("electronNMissingHits",electronNMissingHits,"electronNMissingHits["+TString("nLooseElectron")+"]/I");
    Book("electronTripleCharge",electronTripleCharge,"electronTripleCharge["+TString("nLooseElectron")+"]/I");
    Book("electronCombIso",electronCombIso,"electronCombIso["+TString("nLooseElectron")+"]/F");
    Book("muonSfLoose",muonSfLoose,"muonSfLoose["+TString("nLooseMuon")+"]/F");
    Book("muonSfMedium",muonSfMedium,"muonSfMedium["+TString("nLooseMuon")+"]/F");
    Book("muonSfTight",muonSfTight,"muonSfTight["+TString("nLooseMuon")+"]/F");
    Book("muonSfUnc",muonSfUnc,"muonSfUnc["+TString("nLooseMuon")+"]/F");
    Book("muonSfReco",muonSfReco,"muonSfReco["+TString("nLooseMuon")+"]/F");
    Book("muonD0",muonD0,"muonD0["+TString("nLooseMuon")+"]/F");
    Book("muonDZ",muonDZ,"muonDZ["+TString("nLooseMuon")+"]/F");
    Book("muonIsSoftMuon",muonIsSoftMuon,"muonIsSoftMuon["+TString("nLooseMuon")+"]/I");
    Book("muonCombIso",muonCombIso,"muonCombIso["+TString("nLooseMuon")+"]/F");
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
    Book("genWPlusPt",&genWPlusPt,"genWPlusPt/F");
    Book("genWMinusPt",&genWMinusPt,"genWMinusPt/F");
    Book("genWPlusEta",&genWPlusEta,"genWPlusEta/F");
    Book("genWMinusEta",&genWMinusEta,"genWMinusEta/F");
    Book("looseGenLep1PdgId",&looseGenLep1PdgId,"looseGenLep1PdgId/I");
    Book("looseGenLep2PdgId",&looseGenLep2PdgId,"looseGenLep2PdgId/I");
    Book("looseGenLep3PdgId",&looseGenLep3PdgId,"looseGenLep3PdgId/I");
    Book("looseGenLep4PdgId",&looseGenLep4PdgId,"looseGenLep4PdgId/I");
  }
  if (is_monohiggs||is_hbb) {
    Book("pfmetsig",&pfmetsig,"pfmetsig/F");
    Book("puppimetsig",&puppimetsig,"puppimetsig/F");
    Book("nJotMax",&nJotMax,"nJotMax/I");
    Book("jotM",jotM[0],"jotM["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotM_JESUp",jotM[1],"jotM_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotM_JESDown",jotM[2],"jotM_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotCMVA",jotCMVA[0],"jotCMVA["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotCMVA_JESUp",jotCMVA[1],"jotCMVA_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotCMVA_JESDown",jotCMVA[2],"jotCMVA_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotIso",jotIso[0],"jotIso["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("jotIso_JESUp",jotIso[1],"jotIso_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("jotIso_JESDown",jotIso[2],"jotIso_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("jotQGL",jotQGL[0],"jotQGL["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotQGL_JESUp",jotQGL[1],"jotQGL_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotQGL_JESDown",jotQGL[2],"jotQGL_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotLep1Pt",jotLep1Pt[0],"jotLep1Pt["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotLep1Pt_JESUp",jotLep1Pt[1],"jotLep1Pt_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotLep1Pt_JESDown",jotLep1Pt[2],"jotLep1Pt_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotLep1PtRel",jotLep1PtRel[0],"jotLep1PtRel["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotLep1PtRel_JESUp",jotLep1PtRel[1],"jotLep1PtRel_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotLep1PtRel_JESDown",jotLep1PtRel[2],"jotLep1PtRel_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotLep1DeltaR",jotLep1DeltaR[0],"jotLep1DeltaR["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotLep1DeltaR_JESUp",jotLep1DeltaR[1],"jotLep1DeltaR_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotLep1DeltaR_JESDown",jotLep1DeltaR[2],"jotLep1DeltaR_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotTrk1Pt",jotTrk1Pt[0],"jotTrk1Pt["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotTrk1Pt_JESUp",jotTrk1Pt[1],"jotTrk1Pt_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotTrk1Pt_JESDown",jotTrk1Pt[2],"jotTrk1Pt_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotVtxPt",jotVtxPt[0],"jotVtxPt["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotVtxPt_JESUp",jotVtxPt[1],"jotVtxPt_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotVtxPt_JESDown",jotVtxPt[2],"jotVtxPt_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotVtxMass",jotVtxMass[0],"jotVtxMass["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotVtxMass_JESUp",jotVtxMass[1],"jotVtxMass_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotVtxMass_JESDown",jotVtxMass[2],"jotVtxMass_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotVtx3DVal",jotVtx3DVal[0],"jotVtx3DVal["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotVtx3DVal_JESUp",jotVtx3DVal[1],"jotVtx3DVal_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotVtx3DVal_JESDown",jotVtx3DVal[2],"jotVtx3DVal_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotVtx3DErr",jotVtx3DErr[0],"jotVtx3DErr["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotVtx3DErr_JESUp",jotVtx3DErr[1],"jotVtx3DErr_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotVtx3DErr_JESDown",jotVtx3DErr[2],"jotVtx3DErr_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotVtxNtrk",jotVtxNtrk[0],"jotVtxNtrk["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("jotVtxNtrk_JESUp",jotVtxNtrk[1],"jotVtxNtrk_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("jotVtxNtrk_JESDown",jotVtxNtrk[2],"jotVtxNtrk_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("jotEMF",jotEMF[0],"jotEMF["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotEMF_JESUp",jotEMF[1],"jotEMF_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotEMF_JESDown",jotEMF[2],"jotEMF_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotHF",jotHF[0],"jotHF["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotHF_JESUp",jotHF[1],"jotHF_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotHF_JESDown",jotHF[2],"jotHF_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotNLep",jotNLep[0],"jotNLep["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("jotNLep_JESUp",jotNLep[1],"jotNLep_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("jotNLep_JESDown",jotNLep[2],"jotNLep_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("jotGenPt",jotGenPt[0],"jotGenPt["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotGenPt_JESUp",jotGenPt[1],"jotGenPt_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotGenPt_JESDown",jotGenPt[2],"jotGenPt_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("jotFlav",jotFlav[0],"jotFlav["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("jotFlav_JESUp",jotFlav[1],"jotFlav_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("jotFlav_JESDown",jotFlav[2],"jotFlav_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("hbbjtidx",hbbjtidx[0],"hbbjtidx["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("hbbjtidx_JESUp",hbbjtidx[1],"hbbjtidx_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("hbbjtidx_JESDown",hbbjtidx[2],"hbbjtidx_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/I");
    Book("hbbjtRegFac",hbbjtRegFac[0],"hbbjtRegFac["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("hbbjtRegFac_JESUp",hbbjtRegFac[1],"hbbjtRegFac_JESUp["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("hbbjtRegFac_JESDown",hbbjtRegFac[2],"hbbjtRegFac_JESDown["+TString((is_monohiggs||is_hbb)?"nJotMax":"2")+"]/F");
    Book("fjsjPt",fjsjPt,"fjsjPt["+TString("2")+"]/F");
    Book("fjsjEta",fjsjEta,"fjsjEta["+TString("2")+"]/F");
    Book("fjsjPhi",fjsjPhi,"fjsjPhi["+TString("2")+"]/F");
    Book("fjsjM",fjsjM,"fjsjM["+TString("2")+"]/F");
    Book("fjsjCSV",fjsjCSV,"fjsjCSV["+TString("2")+"]/F");
    Book("fjsjQGL",fjsjQGL,"fjsjQGL["+TString("2")+"]/F");
    Book("hbbpt",&(hbbpt[0]),"hbbpt/F");
    Book("hbbpt_JESUp",&(hbbpt[1]),"hbbpt_JESUp/F");
    Book("hbbpt_JESDown",&(hbbpt[2]),"hbbpt_JESDown/F");
    Book("hbbeta",&(hbbeta[0]),"hbbeta/F");
    Book("hbbeta_JESUp",&(hbbeta[1]),"hbbeta_JESUp/F");
    Book("hbbeta_JESDown",&(hbbeta[2]),"hbbeta_JESDown/F");
    Book("hbbphi",&(hbbphi[0]),"hbbphi/F");
    Book("hbbphi_JESUp",&(hbbphi[1]),"hbbphi_JESUp/F");
    Book("hbbphi_JESDown",&(hbbphi[2]),"hbbphi_JESDown/F");
    Book("hbbm",&(hbbm[0]),"hbbm/F");
    Book("hbbm_JESUp",&(hbbm[1]),"hbbm_JESUp/F");
    Book("hbbm_JESDown",&(hbbm[2]),"hbbm_JESDown/F");
    Book("hbbm_reg",&(hbbm_reg[0]),"hbbm_reg/F");
    Book("hbbm_reg_JESUp",&(hbbm_reg[1]),"hbbm_reg_JESUp/F");
    Book("hbbm_reg_JESDown",&(hbbm_reg[2]),"hbbm_reg_JESDown/F");
    Book("hbbpt_reg",&(hbbpt_reg[0]),"hbbpt_reg/F");
    Book("hbbpt_reg_JESUp",&(hbbpt_reg[1]),"hbbpt_reg_JESUp/F");
    Book("hbbpt_reg_JESDown",&(hbbpt_reg[2]),"hbbpt_reg_JESDown/F");
    Book("nSoft2",&nSoft2,"nSoft2/I");
    Book("nSoft5",&nSoft5,"nSoft5/I");
    Book("nSoft10",&nSoft10,"nSoft10/I");
    Book("hbbCosThetaJJ",&(hbbCosThetaJJ[0]),"hbbCosThetaJJ/F");
    Book("hbbCosThetaJJ_JESUp",&(hbbCosThetaJJ[1]),"hbbCosThetaJJ_JESUp/F");
    Book("hbbCosThetaJJ_JESDown",&(hbbCosThetaJJ[2]),"hbbCosThetaJJ_JESDown/F");
    Book("hbbCosThetaCSJ1",&(hbbCosThetaCSJ1[0]),"hbbCosThetaCSJ1/F");
    Book("hbbCosThetaCSJ1_JESUp",&(hbbCosThetaCSJ1[1]),"hbbCosThetaCSJ1_JESUp/F");
    Book("hbbCosThetaCSJ1_JESDown",&(hbbCosThetaCSJ1[2]),"hbbCosThetaCSJ1_JESDown/F");
    Book("topMassLep1Met",&(topMassLep1Met[0]),"topMassLep1Met/F");
    Book("topMassLep1Met_JESUp",&(topMassLep1Met[1]),"topMassLep1Met_JESUp/F");
    Book("topMassLep1Met_JESDown",&(topMassLep1Met[2]),"topMassLep1Met_JESDown/F");
    Book("topWBosonCosThetaCS",&(topWBosonCosThetaCS[0]),"topWBosonCosThetaCS/F");
    Book("topWBosonCosThetaCS_JESUp",&(topWBosonCosThetaCS[1]),"topWBosonCosThetaCS_JESUp/F");
    Book("topWBosonCosThetaCS_JESDown",&(topWBosonCosThetaCS[2]),"topWBosonCosThetaCS_JESDown/F");
    Book("topWBosonPt",&topWBosonPt,"topWBosonPt/F");
    Book("topWBosonEta",&topWBosonEta,"topWBosonEta/F");
    Book("topWBosonPhi",&topWBosonPhi,"topWBosonPhi/F");
    Book("sumEtSoft1",&sumEtSoft1,"sumEtSoft1/F");
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
    Book("fjGenNumB",&fjGenNumB,"fjGenNumB/I");
  }
  if (is_monotop||is_vbf) {
    Book("jetPt",jetPt,"jetPt["+TString("2")+"]/F");
    Book("jetEta",jetEta,"jetEta["+TString("2")+"]/F");
    Book("jetPhi",jetPhi,"jetPhi["+TString("2")+"]/F");
    Book("jetGenPt",jetGenPt,"jetGenPt["+TString("2")+"]/F");
    Book("jetCSV",jetCSV,"jetCSV["+TString("2")+"]/F");
    Book("jetFlav",jetFlav,"jetFlav["+TString("2")+"]/I");
    Book("jetIsTight",jetIsTight,"jetIsTight["+TString("2")+"]/I");
    Book("jetIsIso",jetIsIso,"jetIsIso["+TString("2")+"]/I");
  }
}