#include "../interface/Module.h"

using namespace pa;
using namespace panda;
using namespace std;
using namespace fastjet; 

ConfigMod::ConfigMod(const Analysis& a_, GeneralTree& gt_, int DEBUG_) :
  BaseModule("config"),
  cfg(a_, DEBUG_),
  analysis(a_),
  gt(gt_),
  bl({"runNumber", "lumiNumber", "eventNumber","weight"})
{

  cfg.isData = analysis.isData;
  utils.eras = new EraHandler(analysis.year);  
  cfg.auxFilePath = analysis.outpath;
  cfg.auxFilePath.ReplaceAll(".root","_aux%i.root");

  // configuration of objects
  if (analysis.ak8)
    cfg.FATJETMATCHDR2 = 0.8*0.8;
  if ((analysis.fatjet || analysis.ak8) || 
      (analysis.recluster || analysis.deep || analysis.deepGen)) {
  }

  if (analysis.recluster || analysis.reclusterGen || 
      analysis.deep || analysis.deepGen || analysis.hbb) {
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    utils.activeArea = new GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
    utils.areaDef = new AreaDefinition(active_area_explicit_ghosts,*(utils.activeArea));
  }

  double radius = 1.5;
  double sdZcut = 0.15;
  double sdBeta = 1.;
  if (analysis.ak8) {
    radius = 0.8;
    sdZcut = 0.1;
    sdBeta = 0.;
  } 
  utils.softDrop = new contrib::SoftDrop(sdBeta,sdZcut,radius);

  if (analysis.deepTracks) {
    cfg.NPFPROPS += 7;
    if (analysis.deepSVs) {
      cfg.NPFPROPS += 3;
    }
  }
  if (analysis.deepExC) {
    cfg.NMAXPF = 250;
    cfg.NGENPROPS = 5;
  }

  if (analysis.hbb) {
    cfg.minJetPt = analysis.ZllHbb ? 20 : 25;
    cfg.minGenFatJetPt = 200;
  }
  if (analysis.vbf || analysis.hbb || analysis.complicatedLeptons) 
    cfg.minBJetPt = 20;
  if (analysis.hbb || analysis.monoh)
    cfg.NJETSAVED = NJET;

  cfg.maxshiftJES = analysis.varyJES ? 
                      (analysis.hbb ? 
                        jes2i(shiftjes::N) : 
                        jes2i(shiftjes::kJESTotalDown) + 1) :
                      1; 

  cfg.ibetas = gt.get_ibetas();
  cfg.Ns = gt.get_Ns();
  cfg.orders = gt.get_orders();

  set_inputBranches();
  set_outputBranches();
}

void ConfigMod::set_inputBranches()
{
  bl.setVerbosity(0);
  TString jetname = analysis.puppiJets ? "puppi" : "chs";

  if (analysis.genOnly) {
    bl += {"genParticles","genReweight","ak4GenJets","genMet","genParticlesU","electrons"};
  } else { 
    bl += {"runNumber", "lumiNumber", "eventNumber", "rho", 
           "isData", "npv", "npvTrue", "weight", "chsAK4Jets", 
           "electrons", "muons", "taus", "photons", 
           "pfMet", "caloMet", "puppiMet", "rawMet", "trkMet", 
           "recoil","metFilters","trkMet"};

    if (analysis.ak8) {
      bl += {jetname+"AK8Jets", "subjets", jetname+"AK8Subjets","Subjets"};
      if (analysis.hbb)
        bl.push_back("ak8GenJets");
    } else if (analysis.fatjet) {
      bl += {jetname+"CA15Jets", "subjets", jetname+"CA15Subjets","Subjets"};
      if (analysis.hbb)
        bl.push_back("ca15GenJets");
    }
    if (analysis.recluster || analysis.bjetRegression || 
        analysis.deep || analysis.hbb || analysis.complicatedPhotons) {
      bl.push_back("pfCandidates");
    }
    if (analysis.deepTracks || analysis.bjetRegression || analysis.hbb) {
      bl += {"tracks","vertices"};
    }

    if (analysis.bjetRegression || analysis.deepSVs)
      bl.push_back("secondaryVertices");

    if (cfg.isData || analysis.mcTriggers) {
      bl.push_back("triggers");
      if (analysis.hbb)
        bl.push_back("triggerObjects");
    }

    if (!cfg.isData) {
      bl += {"genParticles","genReweight","ak4GenJets","genMet"};
      if (analysis.hbb)
        bl.push_back("partons"); 
    }
  }
}

void ConfigMod::set_outputBranches() 
{
  // manipulate the output tree
  if (cfg.isData) {
    std::vector<TString> droppable = {"mcWeight","scale","scaleUp",
                                      "trueGenBosonPt",
                                      "scaleDown","pdf.*","gen.*","sf_.*"};
    gt.RemoveBranches(droppable,{"sf_phoPurity"});
  }
  if (analysis.genOnly) {
    std::vector<TString> keepable = {"mcWeight","scale","scaleUp",
                                     "scaleDown","pdf.*","gen.*",
                                     "sf_tt.*","sf_qcdTT.*",
                                     "trueGenBosonPt","sf_qcd.*","sf_ewk.*",
                                     "nHF.*", "nB.*"};
    if (analysis.deepGen) {
      keepable.push_back("fjECFN.*");
      keepable.push_back("fjTau.*");
    }
    gt.RemoveBranches({".*"},keepable);
  }
  if (!analysis.fatjet && !analysis.ak8) {
    gt.RemoveBranches({"fj.*"});
  }
  if (analysis.complicatedLeptons) {
    gt.RemoveBranches({"genJet.*","puppiU.*","pfU.*","dphipfU.*","dphipuppi.*","jet.*"});
  }
  if (!analysis.hbb) 
    gt.RemoveBranches({".*JES.*"},{".*JESTotal.*"});
  //gt.RemoveBranches({"sf_cmva.*"},{"sf_cmva_Central"});
}

void ConfigMod::readData(TString dirPath)
{
  dirPath += "/";


  // pileup
  utils.openCorr(cNPV,dirPath+"moriond17/normalized_npv.root","data_npv_Wmn",1);
  utils.openCorr(cPU,dirPath+"moriond17/puWeights_80x_37ifb.root","puWeights",1);
  utils.openCorr(cPUUp,dirPath+"moriond17/puWeights_80x_37ifb.root","puWeightsUp",1);
  utils.openCorr(cPUDown,dirPath+"moriond17/puWeights_80x_37ifb.root","puWeightsDown",1);

  if (analysis.complicatedLeptons) {
    // Corrections checked out from Gui's repository on Nov 12, 2017 ~DGH
    // https://github.com/GuillelmoGomezCeballos/MitAnalysisRunII/tree/master/data/80x
    if (analysis.year==2017) {
      utils.openCorr(cEleVeto,
                     dirPath+"leptonic/egammaEffi.txt_EGM2D_runBCDEF_passingVeto94X.root",
                     "EGamma_SF2D",2);
      utils.openCorr(cEleLoose,
                     dirPath+"leptonic/egammaEffi.txt_EGM2D_runBCDEF_passingLoose94X.root",
                     "EGamma_SF2D",2);
      utils.openCorr(cEleMedium,
                     dirPath+"leptonic/egammaEffi.txt_EGM2D_runBCDEF_passingMedium94X.root",
                     "EGamma_SF2D",2);
      utils.openCorr(cEleTight,
                     dirPath+"leptonic/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root",
                     "EGamma_SF2D",2);
      utils.openCorr(cEleMvaWP80,
                     dirPath+"leptonic/egammaEffi.txt_EGM2D_runBCDEF_passingMVA94Xwp80iso.root",
                     "EGamma_SF2D",2);
      utils.openCorr(cEleMvaWP90,
                     dirPath+"leptonic/egammaEffi.txt_EGM2D_runBCDEF_passingMVA94Xwp90iso.root",
                     "EGamma_SF2D",2);
      utils.openCorr(cEleReco,
                     dirPath+"leptonic/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root",
                     "EGamma_SF2D",2);
      // NEED TO UPDATE THESE!
      utils.openCorr(cMuLooseID,
                     dirPath+"leptonic/muon_scalefactors_37ifb.root",
                     "scalefactors_MuonLooseId_Muon",2);
      utils.openCorr(cMuMediumID,
                     dirPath+"leptonic/scalefactors_80x_dylan_37ifb.root",
                     "scalefactors_Medium_Muon",2);
      utils.openCorr(cMuTightID,
                     dirPath+"leptonic/muon_scalefactors_37ifb.root",
                     "scalefactors_TightId_Muon",2);
      utils.openCorr(cMuLooseIso,
                     dirPath+"leptonic/muon_scalefactors_37ifb.root",
                     "scalefactors_Iso_MuonLooseId",2);
      utils.openCorr(cMuMediumIso,
                     dirPath+"leptonic/muon_scalefactors_37ifb.root",
                     "scalefactors_Iso_MuonMediumId",2);
      utils.openCorr(cMuTightIso,
                     dirPath+"leptonic/muon_scalefactors_37ifb.root",
                     "scalefactors_Iso_MuonTightId",2);
      utils.openCorr(cMuReco,
                     dirPath+"leptonic/Tracking_EfficienciesAndSF_BCDEFGH.root",
                     "ratio_eff_eta3_dr030e030_corr",1);
    } else {
      utils.openCorr(cMuLooseID,
                     dirPath+"leptonic/muon_scalefactors_37ifb.root",
                     "scalefactors_MuonLooseId_Muon",2);
      utils.openCorr(cMuMediumID,
                     dirPath+"leptonic/scalefactors_80x_dylan_37ifb.root",
                     "scalefactors_Medium_Muon",2);
      utils.openCorr(cMuTightID,
                     dirPath+"leptonic/muon_scalefactors_37ifb.root",
                     "scalefactors_TightId_Muon",2);
      utils.openCorr(cMuLooseIso,
                     dirPath+"leptonic/muon_scalefactors_37ifb.root",
                     "scalefactors_Iso_MuonLooseId",2);
      utils.openCorr(cMuMediumIso,
                     dirPath+"leptonic/muon_scalefactors_37ifb.root",
                     "scalefactors_Iso_MuonMediumId",2);
      utils.openCorr(cMuTightIso,
                     dirPath+"leptonic/muon_scalefactors_37ifb.root",
                     "scalefactors_Iso_MuonTightId",2);
      utils.openCorr(cMuReco,
                     dirPath+"leptonic/Tracking_EfficienciesAndSF_BCDEFGH.root",
                     "ratio_eff_eta3_dr030e030_corr",1);
      utils.openCorr(cEleVeto,
                     dirPath+"moriond17/scaleFactor_electron_summer16.root",
                     "scaleFactor_electron_vetoid_RooCMSShape_pu_0_100",2);
      utils.openCorr(cEleLoose,
                     dirPath+"leptonic/scalefactors_80x_egpog_37ifb.root",
                     "scalefactors_Loose_Electron",2);
      utils.openCorr(cEleMedium,
                     dirPath+"leptonic/scalefactors_80x_dylan_37ifb.root",
                     "scalefactors_Medium_Electron",2);
      utils.openCorr(cEleTight,
                     dirPath+"leptonic/scalefactors_80x_egpog_37ifb.root",
                     "scalefactors_Tight_Electron",2);
      utils.openCorr(cEleMvaWP90,
                     dirPath+"leptonic/scalefactors_80x_egpog_37ifb.root",
                     "scalefactors_MediumMVA_Electron",2);
      utils.openCorr(cEleMvaWP80,
                     dirPath+"leptonic/scalefactors_80x_egpog_37ifb.root",
                     "scalefactors_TightMVA_Electron",2);
      utils.openCorr(cEleReco,
                     dirPath+"leptonic/scalefactors_80x_egpog_37ifb.root",
                     "scalefactors_Reco_Electron",2);
      utils.openCorr(cEleVeto,
                     dirPath+"moriond17/scaleFactor_electron_summer16.root",
                     "scaleFactor_electron_vetoid_RooCMSShape_pu_0_100",2);
      utils.openCorr(cEleTight,
                     dirPath+"moriond17/scaleFactor_electron_summer16.root",
                     "scaleFactor_electron_tightid_RooCMSShape_pu_0_100",2);
      utils.openCorr(cEleReco,
                     dirPath+"moriond17/scaleFactor_electron_reco_summer16.root",
                     "scaleFactor_electron_reco_RooCMSShape_pu_0_100",2);
      utils.openCorr(cMuLooseID,
                     dirPath+"moriond17/muon_scalefactors_37ifb.root",
                     "scalefactors_MuonLooseId_Muon",2);
      utils.openCorr(cMuLooseIso,
                     dirPath+"moriond17/muon_scalefactors_37ifb.root",
                     "scalefactors_Iso_MuonLooseId",2);
      utils.openCorr(cMuTightID,
                     dirPath+"moriond17/muon_scalefactors_37ifb.root",
                     "scalefactors_TightId_Muon",2);
      utils.openCorr(cMuTightIso,
                     dirPath+"moriond17/muon_scalefactors_37ifb.root",
                     "scalefactors_Iso_MuonTightId",2);
      utils.openCorr(cMuReco,
                     dirPath+"moriond17/Tracking_12p9.root","htrack2",1);
    }
    // EWK corrections 
    utils.openCorr(cWZEwkCorr,
                   dirPath+"leptonic/data.root","hEWKWZCorr",1);
    utils.openCorr(cqqZZQcdCorr,
                   dirPath+"leptonic/data.root","hqqZZKfactor",2);
  } else {
  }
  // Differential Electroweak VH Corrections
  if (analysis.hbb) {
    utils.openCorr(cWmHEwkCorr    ,
                   dirPath+"higgs/Wm_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_rebin"     ,1);
    utils.openCorr(cWmHEwkCorrUp  ,
                   dirPath+"higgs/Wm_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_up_rebin"  ,1);
    utils.openCorr(cWmHEwkCorrDown,
                   dirPath+"higgs/Wm_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_down_rebin",1);
    utils.openCorr(cWpHEwkCorr    ,
                   dirPath+"higgs/Wp_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_rebin"     ,1);
    utils.openCorr(cWpHEwkCorrUp  ,
                   dirPath+"higgs/Wp_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_up_rebin"  ,1);
    utils.openCorr(cWpHEwkCorrDown,
                   dirPath+"higgs/Wp_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_down_rebin",1);
    utils.openCorr(cZnnHEwkCorr    ,
                   dirPath+"higgs/Znn_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_rebin"     ,1);
    utils.openCorr(cZnnHEwkCorrUp  ,
                   dirPath+"higgs/Znn_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_up_rebin"  ,1);
    utils.openCorr(cZnnHEwkCorrDown,
                   dirPath+"higgs/Znn_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_down_rebin",1);
    utils.openCorr(cZllHEwkCorr    ,
                   dirPath+"higgs/Zll_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_rebin"     ,1);
    utils.openCorr(cZllHEwkCorrUp  ,
                   dirPath+"higgs/Zll_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_up_rebin"  ,1);
    utils.openCorr(cZllHEwkCorrDown,
                   dirPath+"higgs/Zll_nloEWK_weight_unnormalized.root",
                   "SignalWeight_nloEWK_down_rebin",1);
  }

  if (analysis.hbb) {
    utils.openCorr(cJetLoosePUID,
                   dirPath+"higgs/puid.root",
                   "puid_76x_loose",2);
  }

  // photons
  utils.openCorr(cPho,
                 dirPath+"moriond17/scalefactors_80x_medium_photon_37ifb.root",
                 "EGamma_SF2D",2);
  utils.openCorr(cPhoFake,
                 dirPath+"moriond17/photonFake.root",
                 "sf",1);

  // triggers
  if (analysis.year==2017) {
    utils.openCorr(cTrigDoubleEleLeg1,
                   dirPath+"leptonic/scalefactors_94x_vhdudes_2017.root",
                   "scalefactors_doubleEleTriggerLeg1",2);
    utils.openCorr(cTrigDoubleEleLeg2,
                   dirPath+"leptonic/scalefactors_94x_vhdudes_2017.root",
                   "scalefactors_doubleEleTriggerLeg2",2);
    utils.openCorr(cTrigMu,
                   dirPath+"trigger_eff/SingleMuonTrigger94x_EfficienciesAndSF_RunBtoF_Nov17Nov2017.root",
                   "IsoMu27_PtEtaBins/abseta_pt_ratio",2);
  } else {
    utils.openCorr(cTrigMu,
                   dirPath+"trigger_eff/muon_trig_Run2016BtoF.root",
                   "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA",2);
    utils.openCorr(cTrigDoubleEleLeg1,
                   dirPath+"trigger_eff/triggers_76x_hww.root",
                   "h2_results_electron_double_leadingleg",2);
    utils.openCorr(cTrigDoubleEleLeg2,
                   dirPath+"trigger_eff/triggers_76x_hww.root",
                   "h2_results_electron_double_trailingleg",2);
    utils.openCorr(cTrigDoubleMuLeg1,
                   dirPath+"trigger_eff/triggers_76x_hww.root",
                   "h2_results_muon_double_leadingleg",2);
    utils.openCorr(cTrigDoubleMuLeg2,
                   dirPath+"trigger_eff/triggers_76x_hww.root",
                   "h2_results_muon_double_trailingleg",2);
  }
  // Currently out of date for 2017
  utils.openCorr(cTrigMET,
                 dirPath+"moriond17/metTriggerEfficiency_recoil_monojet_TH1F.root",
                 "hden_monojet_recoil_clone_passed",1);
  utils.openCorr(cTrigEle,
                 dirPath+"moriond17/eleTrig.root","hEffEtaPt",2);
  utils.openCorr(cTrigPho,
                 dirPath+"moriond17/photonTriggerEfficiency_photon_TH1F.root",
                 "hden_photonpt_clone_passed",1);
  utils.openCorr(cTrigMETZmm,
                 dirPath+"moriond17/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root",
                 "hden_monojet_recoil_clone_passed",1);
  utils.openCorr(cTrigDoubleMuLeg1,
                 dirPath+"trigger_eff/triggers_76x_hww.root",
                 "h2_results_muon_double_leadingleg",2);
  utils.openCorr(cTrigDoubleMuLeg2,
                 dirPath+"trigger_eff/triggers_76x_hww.root",
                 "h2_results_muon_double_trailingleg",2);

  // kfactors
  TFile *fKFactor = analysis.vbf ?
                    new TFile(dirPath+"vbf16/kqcd/kfactor_24bins.root") :
                    new TFile(dirPath+"kfactors.root"); 
  utils.fCorrs[cZNLO] = fKFactor; // just for garbage collection

  TH1F *hZLO    = (TH1F*)fKFactor->Get("ZJets_LO/inv_pt");
  TH1F *hWLO    = (TH1F*)fKFactor->Get("WJets_LO/inv_pt");
  TH1F *hALO    = (TH1F*)fKFactor->Get("GJets_LO/inv_pt_G");

  utils.h1Corrs[cZNLO] = new THCorr1(fKFactor->Get("ZJets_012j_NLO/nominal"));
  utils.h1Corrs[cWNLO] = new THCorr1(fKFactor->Get("WJets_012j_NLO/nominal"));
  utils.h1Corrs[cANLO] = new THCorr1(fKFactor->Get("GJets_1j_NLO/nominal_G"));

  utils.h1Corrs[cZEWK] = new THCorr1(fKFactor->Get("EWKcorr/Z"));
  utils.h1Corrs[cWEWK] = new THCorr1(fKFactor->Get("EWKcorr/W"));
  utils.h1Corrs[cAEWK] = new THCorr1(fKFactor->Get("EWKcorr/photon"));

  utils.h1Corrs[cZEWK]->GetHist()->Divide(utils.h1Corrs[cZNLO]->GetHist());     
  utils.h1Corrs[cWEWK]->GetHist()->Divide(utils.h1Corrs[cWNLO]->GetHist());     
  utils.h1Corrs[cAEWK]->GetHist()->Divide(utils.h1Corrs[cANLO]->GetHist());

  utils.h1Corrs[cZNLO]->GetHist()->Divide(hZLO);    
  utils.h1Corrs[cWNLO]->GetHist()->Divide(hWLO);    
  utils.h1Corrs[cANLO]->GetHist()->Divide(hALO);

  cfg.minGenBosonPt = utils.h1Corrs[cZNLO]->GetHist()->GetBinCenter(1);
  cfg.maxGenBosonPt = utils.h1Corrs[cZNLO]->GetHist()->GetBinCenter(utils.h1Corrs[cZNLO]->GetHist()->GetNbinsX()); 

  utils.openCorr(cANLO2j,dirPath+"moriond17/histo_photons_2jet.root","Func",1);


  if (analysis.vbf) {

    utils.openCorr(cVBF_ZNLO,
                   dirPath+"vbf16/kqcd/mjj/merged_zvv.root","h_kfactors_shape",2);
    utils.openCorr(cVBF_WNLO,
                   dirPath+"vbf16/kqcd/mjj/merged_wlv.root","h_kfactors_shape",2);
    utils.openCorr(cVBF_ZllNLO,
                   dirPath+"vbf16/kqcd/mjj/merged_zll.root","h_kfactors_shape",2);

    utils.openCorr(cVBFTight_ZNLO,
                   dirPath+"vbf16/kqcd/mjj/merged_zvv.root","h_kfactors_cc",1);
    utils.openCorr(cVBFTight_WNLO,
                   dirPath+"vbf16/kqcd/mjj/merged_wlv.root","h_kfactors_cc",1);
    utils.openCorr(cVBFTight_ZllNLO,
                   dirPath+"vbf16/kqcd/mjj/merged_zll.root","h_kfactors_cc",1);


    utils.openCorr(cVBF_EWKZ,
                   dirPath+"vbf16/kewk/kFactor_ZToNuNu_pT_Mjj.root",
                   "TH2F_kFactor",2);
    utils.openCorr(cVBF_EWKW,
                   dirPath+"vbf16/kewk/kFactor_WToLNu_pT_Mjj.root",
                   "TH2F_kFactor",2);

    utils.openCorr(cVBF_TrigMET,
                   dirPath+"vbf16/trig/fit_nmu1.root",
                   "f_eff",3);
    utils.openCorr(cVBF_TrigMETZmm,
                   dirPath+"vbf16/trig/fit_nmu2.root",
                   "f_eff",3);

  }

  utils.openCorr(cBadECALJets,
                 dirPath+"vbf16/hotjets-runBCDEFGH.root",
                 "h2jet",2);


  // btag
  utils.openCorr(cCSVBL,dirPath+"csv/csv_effLoose.root","B",2);
  utils.openCorr(cCSVCL,dirPath+"csv/csv_effLoose.root","C",2);
  utils.openCorr(cCSVLL,dirPath+"csv/csv_effLoose.root","L",2);
  utils.btag = new BTagCorrs(dirPath, analysis, gt); 


  // TODO move these into a getCorr-accessible correction
  // mSD corr
  utils.fMSDcorr = new TFile(dirPath+"/puppiCorr.root");
  utils.puppisd_corrGEN = (TF1*)utils.fMSDcorr->Get("puppiJECcorr_gen");;
  utils.puppisd_corrRECO_cen = (TF1*)utils.fMSDcorr->Get("puppiJECcorr_reco_0eta1v3");
  utils.puppisd_corrRECO_for = (TF1*)utils.fMSDcorr->Get("puppiJECcorr_reco_1v3eta2v5");


  TFile *fcharges = TFile::Open((dirPath + "/deep/charges.root").Data());
  TTree *tcharges = (TTree*)(fcharges->Get("charges"));
  utils.pdgToQ.clear();
  int pdg = 0; float q = 0;
  tcharges->SetBranchAddress("pdgid",&pdg);
  tcharges->SetBranchAddress("q",&q);
  for (unsigned i = 0; i != tcharges->GetEntriesFast(); ++i) {
    tcharges->GetEntry(i);
    utils.pdgToQ[pdg] = q;
  }
  fcharges->Close();
}


void AnalysisMod::initialize(Registry& registry)
{
  if (!on())
    return;
  if (cfg.DEBUG) 
    PDebug("AnalysisMod::initialize", name);
  do_init(registry);
  for (auto* mod: subMods)
    mod->initialize(registry);
}

void AnalysisMod::readData(TString path)
{
  if (!on())
    return;
  if (cfg.DEBUG) 
    PDebug("AnalysisMod::readData", name);
  do_readData(path+"/");
  for (auto* mod: subMods)
    mod->readData(path);
}

// execute DOES NOT cascade down child modules -
// calling subMod execution is left up to the caller
// to allow for more flexibility
void AnalysisMod::execute()
{
  if (!on())
    return;
  if (cfg.DEBUG>1)
    PDebug("AnalysisMod::execute", name);
  do_execute();
  cfg.tr.TriggerEvent("execute "+name);
}

void AnalysisMod::reset()
{
  if (!on())
    return;
  if (cfg.DEBUG) 
    PDebug("AnalysisMod::reset", name);
  do_reset();
  for (auto* mod : subMods)
    mod->reset();
}

void AnalysisMod::terminate()
{
  if (!on())
    return;
  if (cfg.DEBUG) 
    PDebug("AnalysisMod::terminate", name);
  do_terminate();
  for (auto* mod : subMods)
    mod->terminate();
}
