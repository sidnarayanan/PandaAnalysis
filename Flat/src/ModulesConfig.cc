#include "../interface/PandaAnalyzer.h"
#include "TSystem.h"
#include <algorithm>
#include <vector>


// contains any configuration-related method

using namespace panda;
using namespace std;


void PandaAnalyzer::SetOutputFile(TString fOutName) 
{
  fOutPath = fOutName;
  fOut = new TFile(fOutName,"RECREATE");
  fOut->cd();
  tOut = new TTree("events","events");

  fOut->WriteTObject(hDTotalMCWeight);    

  gt->monohiggs      = (analysis->monoh || analysis->hbb);
  gt->vbf            = analysis->vbf;
  gt->fatjet         = (analysis->fatjet || analysis->deepGen);
  gt->leptonic       = analysis->complicatedLeptons;
  gt->photonic       = analysis->complicatedPhotons;
  gt->hfCounting     = analysis->hfCounting;
  gt->btagWeights    = analysis->btagWeights;
  gt->useCMVA        = analysis->useCMVA;

  if (analysis->deep) {
    auxFilePath = fOutName.ReplaceAll(".root","_pf_%u.root");
    IncrementAuxFile();
  }

  if (analysis->deepGen) {
    auxFilePath = fOutName.ReplaceAll(".root","_gen_%u.root");
    IncrementGenAuxFile();
  }

  // fill the signal weights
  for (auto& id : wIDs) 
    gt->signal_weights[id] = 1;

  // Build the input tree here 
  gt->WriteTree(tOut);

  if (DEBUG) PDebug("PandaAnalyzer::SetOutputFile","Created output in "+fOutPath);
}


int PandaAnalyzer::Init(TTree *t, TH1D *hweights, TTree *weightNames)
{
  if (DEBUG) PDebug("PandaAnalyzer::Init","Starting initialization");
  if (!t || !hweights) {
    PError("PandaAnalyzer::Init","Malformed input!");
    return 0;
  }
  tIn = t;

  // create an RNG
  event.rng.setSize(20);

  ////////////////////////////////////////////////////////////////////// 
  // manipulate which branches to read
  event.setStatus(*t, {"!*"}); // turn everything off first

  TString jetname = (analysis->puppiJets) ? "puppi" : "chs";
  panda::utils::BranchList readlist({"runNumber", "lumiNumber", "eventNumber","weight"});
  readlist.setVerbosity(0);

  if (analysis->genOnly) {
    readlist += {"genParticles","genReweight","ak4GenJets","genMet","genParticlesU","electrons"};
  } else { 
    readlist += {"runNumber", "lumiNumber", "eventNumber", "rho", 
                 "isData", "npv", "npvTrue", "weight", "chsAK4Jets", 
                 "electrons", "muons", "taus", "photons", 
                 "pfMet", "caloMet", "puppiMet", "rawMet", "trkMet", 
                 "recoil","metFilters","trkMet"};

    if (analysis->ak8) {
      readlist += {jetname+"AK8Jets", "subjets", jetname+"AK8Subjets","Subjets"};
      if (!isData)
        readlist += {"ak8GenJets"};
    }
    else if (analysis->fatjet) {
      readlist += {jetname+"CA15Jets", "subjets", jetname+"CA15Subjets","Subjets"};
        readlist += {"ca15GenJets"};
    }
    if (analysis->recluster || analysis->bjetRegression || 
        analysis->deep || analysis->hbb || analysis->complicatedPhotons) {
      readlist.push_back("pfCandidates");
    }
    if (analysis->deepTracks || analysis->bjetRegression || analysis->hbb) {
      readlist += {"tracks","vertices"};
    }

    if (analysis->bjetRegression || analysis->deepSVs)
      readlist.push_back("secondaryVertices");

    if (isData || analysis->applyMCTriggers) {
      readlist.push_back("triggers");
    }

    if (!isData) {
      readlist += {"genParticles","genReweight","ak4GenJets","genMet"};
    }
  }

  event.setAddress(*t, readlist); // pass the readlist so only the relevant branches are turned on
  if (DEBUG) PDebug("PandaAnalyzer::Init","Set addresses");

  ////////////////////////////////////////////////////////////////////// 

  // read MC weights
  hDTotalMCWeight = dynamic_cast<TH1D*>(hweights->Clone("hDTotalMCWeight"));
  hDTotalMCWeight->SetDirectory(0);

  if (weightNames && analysis->processType==kSignal) { // hack?
    if (weightNames->GetEntries()!=377 && weightNames->GetEntries()!=22) {
      PError("PandaAnalyzer::Init",
          TString::Format("Reweighting failed because only found %u weights!",
                          unsigned(weightNames->GetEntries())));
      return 1;
    }
    TString *id = new TString();
    weightNames->SetBranchAddress("id",&id);
    unsigned nW = weightNames->GetEntriesFast();
    for (unsigned iW=0; iW!=nW; ++iW) {
      weightNames->GetEntry(iW);
      wIDs.push_back(*id);
    }
  } else if (analysis->processType==kSignal) {
    PError("PandaAnalyzer::Init","This is a signal file, but the weights are missing!");
    return 2;
  }


  ////////////////////////////////////////////////////////////////////// 

  // manipulate the output tree
  gt->RemoveBranches({"ak81.*"}); // unused
  if (isData) {
    std::vector<TString> droppable = {"mcWeight","scale","scaleUp",
                                      "trueGenBosonPt",
                                      "scaleDown","pdf.*","gen.*","sf_.*"};
    gt->RemoveBranches(droppable,{"sf_phoPurity"});
  }
  if (analysis->genOnly) {
    std::vector<TString> keepable = {"mcWeight","scale","scaleUp",
                                     "scaleDown","pdf.*","gen.*",
                                     "sf_tt.*","sf_qcdTT.*",
                                     "trueGenBosonPt","sf_qcd.*","sf_ewk.*",
                                     "nHF.*", "nB.*"};
    if (analysis->deepGen) {
      keepable.push_back("fj1ECFN.*");
      keepable.push_back("fj1Tau.*");
    }
    gt->RemoveBranches({".*"},keepable);
  }
  if (!analysis->fatjet && !analysis->ak8) {
    gt->RemoveBranches({"fj1.*","nAK8jet*"});
  }
  if (analysis->complicatedLeptons) {
    gt->RemoveBranches({"genJet.*","puppiU.*","pfU.*","dphipfU.*","dphipuppi.*","jet1.*","jet2.*"});
  }


  ////////////////////////////////////////////////////////////////////// 

  // give selectors access to the tree
  for (auto* s : selections)
    s->set_gt(gt);

  ////////////////////////////////////////////////////////////////////// 
  
  // remaining configuraiton of objects
  if ((analysis->fatjet || analysis->ak8) || 
      (analysis->recluster || analysis->deep || analysis->deepGen)) {
    double radius = 1.5;
    double sdZcut = 0.15;
    double sdBeta = 1.;
    if (analysis->ak8) {
      radius = 0.8;
      sdZcut = 0.1;
      sdBeta = 0.;
      jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm,radius);
      jetDefKt = new fastjet::JetDefinition(fastjet::kt_algorithm,radius);
    } else {
      radius = 1.5;
      sdZcut = 0.15;
      sdBeta = 1.;
      jetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm,radius);
      jetDefKt = new fastjet::JetDefinition(fastjet::kt_algorithm,radius);
    }
    softDrop = new fastjet::contrib::SoftDrop(sdBeta,sdZcut,radius);
    tauN = new fastjet::contrib::Njettiness(fastjet::contrib::OnePass_KT_Axes(), 
                                            fastjet::contrib::NormalizedMeasure(1., radius));

    if (analysis->deepGen) {
      ecfcalc = new pandaecf::Calculator();
      if (analysis->deepGenGrid) {
        // grid = new ParticleGridder(250,157,5); // 0.02x0.02
        // grid = new ParticleGridder(1000,628,5); // 0.005x0.005
        grid = new ParticleGridder(2500,1570,5); // 0.002x0.002
        // grid->_on = false;
        grid->_etaphi = false;
      }
    }
  }

  if (analysis->recluster || analysis->reclusterGen || analysis->deep || analysis->deepGen || analysis->hbb) {
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    activeArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
    areaDef = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*activeArea);
  }

  if (analysis->deepTracks) {
    NPFPROPS += 7;
    if (analysis->deepSVs) {
      NPFPROPS += 3;
    }
  }
  if (analysis->deepExC) {
    NMAXPF = 250;
    NGENPROPS = 5;
  }

  if (analysis->reclusterGen || analysis->deepExC) {
    double radius = 0.4;
    jetDefGen = new fastjet::JetDefinition(fastjet::antikt_algorithm,radius);
  }
  if (analysis->hbb)
    softTrackJetDefinition = new fastjet::JetDefinition(fastjet::antikt_algorithm,0.4);

  // Custom jet pt threshold
  if (analysis->hbb) {
    jetPtThreshold = analysis->ZllHbb? 20:25;
    genFatJetMinPt = 200;
  }
  if (analysis->vbf || analysis->hbb || analysis->complicatedLeptons) 
    bJetPtThreshold = 20;

  if (DEBUG) PDebug("PandaAnalyzer::Init","Finished configuration");

  return 0;
}


void PandaAnalyzer::OpenCorrection(CorrectionType ct, TString fpath, TString hname, int dim) 
{
  fCorrs[ct] = TFile::Open(fpath);
  if (dim==1) 
    h1Corrs[ct] = new THCorr1(fCorrs[ct]->Get(hname));
  else if (dim==2)
    h2Corrs[ct] = new THCorr2(fCorrs[ct]->Get(hname));
  else 
    f1Corrs[ct] = new TF1Corr((TF1*)fCorrs[ct]->Get(hname));
}


void PandaAnalyzer::SetDataDir(const char *s) 
{
  TString dirPath(s);
  dirPath += "/";

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Starting loading of data");

  // pileup
  OpenCorrection(cNPV,dirPath+"moriond17/normalized_npv.root","data_npv_Wmn",1);
  OpenCorrection(cPU,dirPath+"moriond17/puWeights_80x_37ifb.root","puWeights",1);
  OpenCorrection(cPUUp,dirPath+"moriond17/puWeights_80x_37ifb.root","puWeightsUp",1);
  OpenCorrection(cPUDown,dirPath+"moriond17/puWeights_80x_37ifb.root","puWeightsDown",1);

  if (analysis->complicatedLeptons) {
    // Corrections checked out from Gui's repository on Nov 12, 2017 ~DGH
    // https://github.com/GuillelmoGomezCeballos/MitAnalysisRunII/tree/master/data/80x
    OpenCorrection(cMuLooseID,
                   dirPath+"leptonic/muon_scalefactors_37ifb.root",
                   "scalefactors_MuonLooseId_Muon",2);
    OpenCorrection(cMuMediumID,
                   dirPath+"leptonic/scalefactors_80x_dylan_37ifb.root",
                   "scalefactors_Medium_Muon",2);
    OpenCorrection(cMuTightID,
                   dirPath+"leptonic/muon_scalefactors_37ifb.root",
                   "scalefactors_TightId_Muon",2);
    OpenCorrection(cMuLooseIso,
                   dirPath+"leptonic/muon_scalefactors_37ifb.root",
                   "scalefactors_Iso_MuonLooseId",2);
    OpenCorrection(cMuMediumIso,
                   dirPath+"leptonic/muon_scalefactors_37ifb.root",
                   "scalefactors_Iso_MuonMediumId",2);
    OpenCorrection(cMuTightIso,
                   dirPath+"leptonic/muon_scalefactors_37ifb.root",
                   "scalefactors_Iso_MuonTightId",2);
    OpenCorrection(cMuReco,
                   dirPath+"leptonic/Tracking_EfficienciesAndSF_BCDEFGH.root",
                   "ratio_eff_eta3_dr030e030_corr",1);
    OpenCorrection(cEleVeto,
                   dirPath+"moriond17/scaleFactor_electron_summer16.root",
                   "scaleFactor_electron_vetoid_RooCMSShape_pu_0_100",2);
    OpenCorrection(cEleLoose,
                   dirPath+"leptonic/scalefactors_80x_egpog_37ifb.root",
                   "scalefactors_Loose_Electron",2);
    OpenCorrection(cEleMedium,
                   dirPath+"leptonic/scalefactors_80x_dylan_37ifb.root",
                   "scalefactors_Medium_Electron",2);
    OpenCorrection(cEleTight,
                   dirPath+"leptonic/scalefactors_80x_egpog_37ifb.root",
                   "scalefactors_Tight_Electron",2);
    OpenCorrection(cEleMvaWP90,
                   dirPath+"leptonic/scalefactors_80x_egpog_37ifb.root",
                   "scalefactors_MediumMVA_Electron",2);
    OpenCorrection(cEleMvaWP80,
                   dirPath+"leptonic/scalefactors_80x_egpog_37ifb.root",
                   "scalefactors_TightMVA_Electron",2);
    OpenCorrection(cEleReco,
                   dirPath+"leptonic/scalefactors_80x_egpog_37ifb.root",
                   "scalefactors_Reco_Electron",2);
    // EWK corrections 
    OpenCorrection(cWZEwkCorr,
                   dirPath+"leptonic/data.root","hEWKWZCorr",1);
    OpenCorrection(cqqZZQcdCorr,
                   dirPath+"leptonic/data.root","hqqZZKfactor",2);

    // TO DO: Hard coded to 2016 rochester corrections for now, need to do this in a better way later
    rochesterCorrection = new RoccoR(Form("%s/rcdata.2016.v3",dirPath.Data()));
  } else {
    OpenCorrection(cEleVeto,
                   dirPath+"moriond17/scaleFactor_electron_summer16.root",
                   "scaleFactor_electron_vetoid_RooCMSShape_pu_0_100",2);
    OpenCorrection(cEleTight,
                   dirPath+"moriond17/scaleFactor_electron_summer16.root",
                   "scaleFactor_electron_tightid_RooCMSShape_pu_0_100",2);
    OpenCorrection(cEleReco,
                   dirPath+"moriond17/scaleFactor_electron_reco_summer16.root",
                   "scaleFactor_electron_reco_RooCMSShape_pu_0_100",2);
    OpenCorrection(cMuLooseID,
                   dirPath+"moriond17/muon_scalefactors_37ifb.root",
                   "scalefactors_MuonLooseId_Muon",2);
    OpenCorrection(cMuLooseIso,
                   dirPath+"moriond17/muon_scalefactors_37ifb.root",
                   "scalefactors_Iso_MuonLooseId",2);
    OpenCorrection(cMuTightID,
                   dirPath+"moriond17/muon_scalefactors_37ifb.root",
                   "scalefactors_TightId_Muon",2);
    OpenCorrection(cMuTightIso,
                   dirPath+"moriond17/muon_scalefactors_37ifb.root",
                   "scalefactors_Iso_MuonTightId",2);
    OpenCorrection(cMuReco,
                   dirPath+"moriond17/Tracking_12p9.root","htrack2",1);
  }
  // Differential Electroweak VH Corrections
  if (analysis->hbb) {
    OpenCorrection(cWmHEwkCorr    ,
                   dirPath+"higgs/Wm_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_rebin"     ,1);
    OpenCorrection(cWmHEwkCorrUp  ,
                   dirPath+"higgs/Wm_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_up_rebin"  ,1);
    OpenCorrection(cWmHEwkCorrDown,
                   dirPath+"higgs/Wm_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_down_rebin",1);
    OpenCorrection(cWpHEwkCorr    ,
                   dirPath+"higgs/Wp_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_rebin"     ,1);
    OpenCorrection(cWpHEwkCorrUp  ,
                   dirPath+"higgs/Wp_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_up_rebin"  ,1);
    OpenCorrection(cWpHEwkCorrDown,
                   dirPath+"higgs/Wp_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_down_rebin",1);
    OpenCorrection(cZnnHEwkCorr    ,
                   dirPath+"higgs/Znn_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_rebin"     ,1);
    OpenCorrection(cZnnHEwkCorrUp  ,
                   dirPath+"higgs/Znn_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_up_rebin"  ,1);
    OpenCorrection(cZnnHEwkCorrDown,
                   dirPath+"higgs/Znn_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_down_rebin",1);
    OpenCorrection(cZllHEwkCorr    ,
                   dirPath+"higgs/Zll_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_rebin"     ,1);
    OpenCorrection(cZllHEwkCorrUp  ,
                   dirPath+"higgs/Zll_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_up_rebin"  ,1);
    OpenCorrection(cZllHEwkCorrDown,
                   dirPath+"higgs/Zll_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_down_rebin",1);
  }

  // photons
  OpenCorrection(cPho,
                 dirPath+"moriond17/scalefactors_80x_medium_photon_37ifb.root",
                 "EGamma_SF2D",2);

  // triggers
  OpenCorrection(cTrigMET,
                 dirPath+"moriond17/metTriggerEfficiency_recoil_monojet_TH1F.root",
                 "hden_monojet_recoil_clone_passed",1);
  OpenCorrection(cTrigEle,
                 dirPath+"moriond17/eleTrig.root","hEffEtaPt",2);
  OpenCorrection(cTrigMu,
                 dirPath+"trigger_eff/muon_trig_Run2016BtoF.root",
                 "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA",2);
  OpenCorrection(cTrigPho,
                 dirPath+"moriond17/photonTriggerEfficiency_photon_TH1F.root",
                 "hden_photonpt_clone_passed",1);
  OpenCorrection(cTrigMETZmm,
                 dirPath+"moriond17/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root",
                 "hden_monojet_recoil_clone_passed",1);

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded scale factors");

  // kfactors
  TFile *fKFactor = 0;
  if (analysis->vbf)
    fKFactor = new TFile(dirPath+"vbf16/kqcd/kfactor_24bins.root"); 
  else
    fKFactor = new TFile(dirPath+"kfactors.root"); 
  fCorrs[cZNLO] = fKFactor; // just for garbage collection

  TH1F *hZLO    = (TH1F*)fKFactor->Get("ZJets_LO/inv_pt");
  TH1F *hWLO    = (TH1F*)fKFactor->Get("WJets_LO/inv_pt");
  TH1F *hALO    = (TH1F*)fKFactor->Get("GJets_LO/inv_pt_G");

  h1Corrs[cZNLO] = new THCorr1(fKFactor->Get("ZJets_012j_NLO/nominal"));
  h1Corrs[cWNLO] = new THCorr1(fKFactor->Get("WJets_012j_NLO/nominal"));
  h1Corrs[cANLO] = new THCorr1(fKFactor->Get("GJets_1j_NLO/nominal_G"));

  h1Corrs[cZEWK] = new THCorr1(fKFactor->Get("EWKcorr/Z"));
  h1Corrs[cWEWK] = new THCorr1(fKFactor->Get("EWKcorr/W"));
  h1Corrs[cAEWK] = new THCorr1(fKFactor->Get("EWKcorr/photon"));

  h1Corrs[cZEWK]->GetHist()->Divide(h1Corrs[cZNLO]->GetHist());     
  h1Corrs[cWEWK]->GetHist()->Divide(h1Corrs[cWNLO]->GetHist());     
  h1Corrs[cAEWK]->GetHist()->Divide(h1Corrs[cANLO]->GetHist());

  h1Corrs[cZNLO]->GetHist()->Divide(hZLO);    
  h1Corrs[cWNLO]->GetHist()->Divide(hWLO);    
  h1Corrs[cANLO]->GetHist()->Divide(hALO);

  OpenCorrection(cANLO2j,dirPath+"moriond17/histo_photons_2jet.root","Func",1);

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded k factors");

  if (analysis->vbf) {

    OpenCorrection(cVBF_ZNLO,
                   dirPath+"vbf16/kqcd/mjj/merged_zvv.root","h_kfactors_shape",2);
    OpenCorrection(cVBF_WNLO,
                   dirPath+"vbf16/kqcd/mjj/merged_wlv.root","h_kfactors_shape",2);
    OpenCorrection(cVBF_ZllNLO,
                   dirPath+"vbf16/kqcd/mjj/merged_zll.root","h_kfactors_shape",2);

    OpenCorrection(cVBFTight_ZNLO,
                   dirPath+"vbf16/kqcd/mjj/merged_zvv.root","h_kfactors_cc",1);
    OpenCorrection(cVBFTight_WNLO,
                   dirPath+"vbf16/kqcd/mjj/merged_wlv.root","h_kfactors_cc",1);
    OpenCorrection(cVBFTight_ZllNLO,
                   dirPath+"vbf16/kqcd/mjj/merged_zll.root","h_kfactors_cc",1);

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded VBF k factors");

    OpenCorrection(cVBF_EWKZ,
                   dirPath+"vbf16/kewk/kFactor_ZToNuNu_pT_Mjj.root",
                   "TH2F_kFactor",2);
    OpenCorrection(cVBF_EWKW,
                   dirPath+"vbf16/kewk/kFactor_WToLNu_pT_Mjj.root",
                   "TH2F_kFactor",2);

    OpenCorrection(cVBF_TrigMET,
                   dirPath+"vbf16/trig/fit_nmu1.root",
                   "f_eff",3);
    OpenCorrection(cVBF_TrigMETZmm,
                   dirPath+"vbf16/trig/fit_nmu2.root",
                   "f_eff",3);

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded VBF k factors");
  }

  OpenCorrection(cBadECALJets,
                 dirPath+"vbf16/hotjets-runBCDEFGH.root",
                 "h2jet",2);


  if (analysis->btagSFs) {
    // btag SFs
    btagCalib = new BTagCalibration("csvv2",(dirPath+"moriond17/CSVv2_Moriond17_B_H.csv").Data());
    btagReaders[bJetL] = new BTagCalibrationReader(BTagEntry::OP_LOOSE,"central",{"up","down"});
    btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_B,"comb");
    btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_C,"comb");
    btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_UDSG,"incl");

    sj_btagCalib = new BTagCalibration("csvv2",(dirPath+"moriond17/subjet_CSVv2_Moriond17_B_H.csv").Data());
    btagReaders[bSubJetL] = new BTagCalibrationReader(BTagEntry::OP_LOOSE,"central",{"up","down"});
    btagReaders[bSubJetL]->load(*sj_btagCalib,BTagEntry::FLAV_B,"lt");
    btagReaders[bSubJetL]->load(*sj_btagCalib,BTagEntry::FLAV_C,"lt");
    btagReaders[bSubJetL]->load(*sj_btagCalib,BTagEntry::FLAV_UDSG,"incl");

    btagReaders[bJetM] = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"});
    btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_B,"comb");
    btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_C,"comb");
    btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_UDSG,"incl");
    
    btagReaders[bSubJetM] = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"});
    btagReaders[bSubJetM]->load(*sj_btagCalib,BTagEntry::FLAV_B,"lt");
    btagReaders[bSubJetM]->load(*sj_btagCalib,BTagEntry::FLAV_C,"lt");
    btagReaders[bSubJetM]->load(*sj_btagCalib,BTagEntry::FLAV_UDSG,"incl");

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded btag SFs");
  } 
  if (analysis->btagWeights) {
    if (analysis->useCMVA) 
      cmvaReweighter = new CSVHelper(
            "PandaAnalysis/data/csvweights/cmva_rwt_fit_hf_v0_final_2017_3_29.root", 
            "PandaAnalysis/data/csvweights/cmva_rwt_fit_lf_v0_final_2017_3_29.root", 
            5
          );
    else
      csvReweighter  = new CSVHelper(
            "PandaAnalysis/data/csvweights/csv_rwt_fit_hf_v2_final_2017_3_29test.root", 
            "PandaAnalysis/data/csvweights/csv_rwt_fit_lf_v2_final_2017_3_29test.root", 
            5
          );
  }

  // bjet regression
  if (analysis->bjetRegression) {
    bjetreg_vars = new float[10];
    bjetregReader = new TMVA::Reader("!Color:!Silent");

    bjetregReader->AddVariable("jetPt[hbbjtidx[0]]",&bjetreg_vars[0]);
    bjetregReader->AddVariable("nJot",&bjetreg_vars[1]);
    bjetregReader->AddVariable("jetEta[hbbjtidx[0]]",&bjetreg_vars[2]);
    bjetregReader->AddVariable("jetE[hbbjtidx[0]]",&bjetreg_vars[3]);
    bjetregReader->AddVariable("npv",&bjetreg_vars[4]);
    bjetregReader->AddVariable("jetLeadingTrkPt[hbbjtidx[0]]",&bjetreg_vars[5]);
    bjetregReader->AddVariable("jetLeadingLepPt[hbbjtidx[0]]",&bjetreg_vars[6]);
    bjetregReader->AddVariable("jetNLep[hbbjtidx[0]]",&bjetreg_vars[7]);
    bjetregReader->AddVariable("jetEMFrac[hbbjtidx[0]]",&bjetreg_vars[8]);
    bjetregReader->AddVariable("jetHadFrac[hbbjtidx[0]]",&bjetreg_vars[9]);

    gSystem->Exec(
        Form("wget -q -O %s/trainings/bjet_regression_v0.weights.xml http://t3serv001.mit.edu/~snarayan/pandadata/trainings/bjet_regression_v0.weights.xml",dirPath.Data())
      );
    bjetregReader->BookMVA( "BDT method", dirPath+"trainings/bjet_regression_v0.weights.xml" );

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded bjet regression weights");
  }


  if (analysis->monoh || analysis->hbb) {
    // mSD corr
    MSDcorr = new TFile(dirPath+"/puppiCorr.root");
    puppisd_corrGEN = (TF1*)MSDcorr->Get("puppiJECcorr_gen");;
    puppisd_corrRECO_cen = (TF1*)MSDcorr->Get("puppiJECcorr_reco_0eta1v3");
    puppisd_corrRECO_for = (TF1*)MSDcorr->Get("puppiJECcorr_reco_1v3eta2v5");

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded mSD correction");
  }

  if (analysis->rerunJES) {
    TString jecV = "V4", jecReco = "23Sep2016"; 
    TString jecVFull = jecReco+jecV;
    ak8UncReader["MC"] = new JetCorrectionUncertainty(
         (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_Uncertainty_AK8PFPuppi.txt").Data()
      );
    std::vector<TString> eraGroups = {"BCD","EF","G","H"};
    for (auto e : eraGroups) {
      ak8UncReader["data"+e] = new JetCorrectionUncertainty(
           (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_Uncertainty_AK8PFPuppi.txt").Data()
        );
    }

    ak8JERReader = new JERReader(dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_SF_AK8PFPuppi.txt",
                                 dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt");


    ak4UncReader["MC"] = new JetCorrectionUncertainty(
         (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_Uncertainty_AK4PFPuppi.txt").Data()
      );
    for (auto e : eraGroups) {
      ak4UncReader["data"+e] = new JetCorrectionUncertainty(
           (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_Uncertainty_AK4PFPuppi.txt").Data()
        );
    }

    ak4JERReader = new JERReader(dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_SF_AK4PFPuppi.txt",
                                 dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_PtResolution_AK4PFPuppi.txt");

    std::vector<JetCorrectorParameters> params = {
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L1FastJet_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L2Relative_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L3Absolute_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L2L3Residual_AK4PFPuppi.txt").Data())
    };
    ak4ScaleReader["MC"] = new FactorizedJetCorrector(params);
    if (DEBUG>1) PDebug("PandaAnalyzer::SetDataDir","Loaded JES for AK4 MC");
    for (auto e : eraGroups) {
      params = {
        JetCorrectorParameters(
          (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L1FastJet_AK4PFPuppi.txt").Data()),
        JetCorrectorParameters(
          (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L2Relative_AK4PFPuppi.txt").Data()),
        JetCorrectorParameters(
          (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L3Absolute_AK4PFPuppi.txt").Data()),
        JetCorrectorParameters(
          (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L2L3Residual_AK4PFPuppi.txt").Data())
      };
      ak4ScaleReader["data"+e] = new FactorizedJetCorrector(params);
      if (DEBUG>1) PDebug("PandaAnalyzer::SetDataDir","Loaded JES for AK4 "+e);
    }

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded JES/R");
  }


  if (analysis->deepGen) {
    TFile *fcharges = TFile::Open((dirPath + "/deep/charges.root").Data());
    TTree *tcharges = (TTree*)(fcharges->Get("charges"));
    pdgToQ.clear();
    int pdg = 0; float q = 0;
    tcharges->SetBranchAddress("pdgid",&pdg);
    tcharges->SetBranchAddress("q",&q);
    for (unsigned i = 0; i != tcharges->GetEntriesFast(); ++i) {
      tcharges->GetEntry(i);
      pdgToQ[pdg] = q;
    }
    fcharges->Close();
  }

}


void PandaAnalyzer::AddGoodLumiRange(int run, int l0, int l1) 
{
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) { // don't know about this run yet
    std::vector<LumiRange> newLumiList;
    newLumiList.emplace_back(l0,l1);
    goodLumis[run] = newLumiList;
  } else {
    run_->second.emplace_back(l0,l1);
  }
}


bool PandaAnalyzer::PassGoodLumis(int run, int lumi) 
{
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) {
    // matched no run
    if (DEBUG) 
      PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i",run));
    return false;
  }

  // found the run, now look for a lumi range
  for (auto &range : run_->second) {
    if (range.Contains(lumi)) {
      if (DEBUG) 
        PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Accepting run=%i, lumi=%i",run,lumi));
      return true;
    }
  }

  // matched no lumi range
  if (DEBUG) 
    PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i, lumi=%i",run,lumi));
  return false;
}





