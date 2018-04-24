#include "../interface/Module.h"

using namespace pa;
using namespace panda;
using namespace std;

ConfigMod::ConfigMod(const Analysis& a_, const GeneralTree& gt, int DEBUG_) :
  analysis(a_),
  cfg(analysis, DEBUG_)
{

  cfg.isData = analysis.isData;

  // remaining configuraiton of objects
  if ((analysis.fatjet || analysis.ak8) || 
      (analysis.recluster || analysis.deep || analysis.deepGen)) {
    double radius = 1.5;
    double sdZcut = 0.15;
    double sdBeta = 1.;
    if (analysis.ak8) {
      radius = 0.8;
      sdZcut = 0.1;
      sdBeta = 0.;
      utils.jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm,radius);
      utils.jetDefKt = new fastjet::JetDefinition(fastjet::kt_algorithm,radius);
    } else {
      radius = 1.5;
      sdZcut = 0.15;
      sdBeta = 1.;
      utils.jetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm,radius);
      utils.jetDefKt = new fastjet::JetDefinition(fastjet::kt_algorithm,radius);
    }
    utils.softDrop = new fastjet::contrib::SoftDrop(sdBeta,sdZcut,radius);
    utils.tauN = new fastjet::contrib::Njettiness(fastjet::contrib::OnePass_KT_Axes(), 
                                                   fastjet::contrib::NormalizedMeasure(1., radius));

    if (analysis.deepGen) {
      utils.ecfcalc = new ECFCalculator();
      if (analysis.deepGenGrid) {
        utils.grid = new ParticleGridder(2500,1570,5); // 0.002x0.002
        utils.grid->_etaphi = false;
      }
    }
  }

  if (analysis.recluster || analysis.reclusterGen || analysis.deep || analysis.deepGen || analysis.hbb) {
    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    utils.activeArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
    utils.areaDef = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*(utils.activeArea));
  }

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

  if (analysis.reclusterGen || analysis.deepExC) {
    double radius = 0.4;
    utils.jetDefGen = new fastjet::JetDefinition(fastjet::antikt_algorithm,radius);
  }
  if (analysis.hbb)
    utils.softTrackJetDefinition = new fastjet::JetDefinition(fastjet::antikt_algorithm,0.4);

  // Custom jet pt threshold
  if (analysis.hbb) {
    cfg.minJetPt = analysis.ZllHbb ? 20 : 25;
    cfg.minGenFatJetPt = 200;
  }
  if (analysis.vbf || analysis.hbb || analysis.complicatedLeptons) 
    cfg.minBJetPt = 20;

  cfg.maxshiftJES = analysis.varyJES ? jes2i(shiftjes::N) : 1; 

  cfg.ibetas = gt->get_ibetas();
  cfg.Ns = gt->get_Ns();
  cfg.orders = gt->get_orders();
}


void ConfigMod::SetDataDir(TString dirPath)
{
  dirPath += "/";

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Starting loading of data");

  // pileup
  utils.openCorr(cNPV,dirPath+"moriond17/normalized_npv.root","data_npv_Wmn",1);
  utils.openCorr(cPU,dirPath+"moriond17/puWeights_80x_37ifb.root","puWeights",1);
  utils.openCorr(cPUUp,dirPath+"moriond17/puWeights_80x_37ifb.root","puWeightsUp",1);
  utils.openCorr(cPUDown,dirPath+"moriond17/puWeights_80x_37ifb.root","puWeightsDown",1);

  if (analysis.complicatedLeptons) {
    // Corrections checked out from Gui's repository on Nov 12, 2017 ~DGH
    // https://github.com/GuillelmoGomezCeballos/MitAnalysisRunII/tree/master/data/80x
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
    // EWK corrections 
    utils.openCorr(cWZEwkCorr,
                   dirPath+"leptonic/data.root","hEWKWZCorr",1);
    utils.openCorr(cqqZZQcdCorr,
                   dirPath+"leptonic/data.root","hqqZZKfactor",2);

    // TO DO: Hard coded to 2016 rochester corrections for now, need to do this in a better way later
    utils.rochesterCorrection = new RoccoR(Form("%s/rcdata.2016.v3",dirPath.Data()));
  } else {
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

  // triggers
  utils.openCorr(cTrigMET,
                 dirPath+"moriond17/metTriggerEfficiency_recoil_monojet_TH1F.root",
                 "hden_monojet_recoil_clone_passed",1);
  utils.openCorr(cTrigEle,
                 dirPath+"moriond17/eleTrig.root","hEffEtaPt",2);
  utils.openCorr(cTrigMu,
                 dirPath+"trigger_eff/muon_trig_Run2016BtoF.root",
                 "IsoMu24_OR_IsoTkMu24_PtEtaBins/efficienciesDATA/abseta_pt_DATA",2);
  utils.openCorr(cTrigPho,
                 dirPath+"moriond17/photonTriggerEfficiency_photon_TH1F.root",
                 "hden_photonpt_clone_passed",1);
  utils.openCorr(cTrigMETZmm,
                 dirPath+"moriond17/metTriggerEfficiency_zmm_recoil_monojet_TH1F.root",
                 "hden_monojet_recoil_clone_passed",1);

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded scale factors");

  // kfactors
  TFile *fKFactor = 0;
  if (analysis.vbf)
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

  utils.openCorr(cANLO2j,dirPath+"moriond17/histo_photons_2jet.root","Func",1);

  if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded k factors");

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

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded VBF k factors");

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

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded VBF k factors");
  }

  utils.openCorr(cBadECALJets,
                 dirPath+"vbf16/hotjets-runBCDEFGH.root",
                 "h2jet",2);


  if (analysis.btagSFs) {
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
  if (analysis.btagWeights) {
    if (analysis.useCMVA) 
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
  if (analysis.bjetRegression) {
    bjetreg_vars = new float[15];
    bjetregReader = new TMVA::Reader("!Color:!Silent");
    bjetregReader->AddVariable("jetPt[hbbjtidx[0]]",
                               &bjetreg_vars[ 0] );
    bjetregReader->AddVariable("jetEta[hbbjtidx[0]]",
                               &bjetreg_vars[ 1] );
    bjetregReader->AddVariable("jetLeadingTrkPt[hbbjtidx[0]]",
                               &bjetreg_vars[ 2] );
    bjetregReader->AddVariable("jetLeadingLepPt[hbbjtidx[0]]",
                               &bjetreg_vars[ 3] );
    bjetregReader->AddVariable("jetEMFrac[hbbjtidx[0]]",
                               &bjetreg_vars[ 4] );
    bjetregReader->AddVariable("jetHadFrac[hbbjtidx[0]]",
                               &bjetreg_vars[ 5] );
    bjetregReader->AddVariable("jetLeadingLepDeltaR[hbbjtidx[0]]",
                               &bjetreg_vars[ 6] );
    bjetregReader->AddVariable("jetLeadingLepPtRel[hbbjtidx[0]]",
                               &bjetreg_vars[ 7] );
    bjetregReader->AddVariable("jetvtxPt[hbbjtidx[0]]",
                               &bjetreg_vars[ 8] );
    bjetregReader->AddVariable("jetvtxMass[hbbjtidx[0]]",
                               &bjetreg_vars[ 9] );
    bjetregReader->AddVariable("jetvtx3Dval[hbbjtidx[0]]",
                               &bjetreg_vars[10] );
    bjetregReader->AddVariable("jetvtx3Derr[hbbjtidx[0]]",
                               &bjetreg_vars[11] );
    bjetregReader->AddVariable("jetvtxNtrk[hbbjtidx[0]]",
                               &bjetreg_vars[12] );
    bjetregReader->AddVariable("evalEt(jetPt[hbbjtidx[0]],jetEta[hbbjtidx[0]],jetPhi[hbbjtidx[0]],jetE[hbbjtidx[0]])",
                               &bjetreg_vars[13] );
    bjetregReader->AddVariable("evalMt(jetPt[hbbjtidx[0]],jetEta[hbbjtidx[0]],jetPhi[hbbjtidx[0]],jetE[hbbjtidx[0]])",
                               &bjetreg_vars[14] );

    gSystem->Exec(
        Form("wget -nv -O %s/trainings/bjet_regression_v1_fromBenedikt.weights.xml http://t3serv001.mit.edu/~dhsu/pandadata/trainings/bjet_regression_v1_fromBenedikt.weights.xml",dirPath.Data())
      );
    bjetregReader->BookMVA( "BDT method", dirPath+"trainings/bjet_regression_v1_fromBenedikt.weights.xml" );

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded bjet regression weights");
  }


  if (analysis.monoh || analysis.hbb) {
    // mSD corr
    MSDcorr = new TFile(dirPath+"/puppiCorr.root");
    puppisd_corrGEN = (TF1*)MSDcorr->Get("puppiJECcorr_gen");;
    puppisd_corrRECO_cen = (TF1*)MSDcorr->Get("puppiJECcorr_reco_0eta1v3");
    puppisd_corrRECO_for = (TF1*)MSDcorr->Get("puppiJECcorr_reco_1v3eta2v5");

    if (DEBUG) PDebug("PandaAnalyzer::SetDataDir","Loaded mSD correction");
  }

  if (analysis.rerunJES) {
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


  if (analysis.deepGen) {
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
