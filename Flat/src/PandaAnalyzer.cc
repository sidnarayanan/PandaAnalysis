#include "../interface/PandaAnalyzer.h"
#include "TMath.h"
#include <algorithm>
#include <vector>


using namespace panda;
using namespace std;
using namespace pa;

#define ADDMOD(X) mods_all.back()->addSubMod<X>()


PandaAnalyzer::PandaAnalyzer(Analysis* a, int debug_/*=0*/) :
  Analyzer("PandaAnalyzer", a, debug_),
  cfgmod(analysis, gt, DEBUG),
  wIDs(v_make_shared<TString>())
{
  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Calling constructor");

  Config& cfg = cfgmod.cfg;
  Utils& utils = cfgmod.utils;

  gblmod = new GlobalMod(event, cfg, utils, gt);
  mods_all.emplace_back(gblmod);

  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Adding AnalysisMods");

  // Define analyses
  preselmod = new ContainerMod("pre-sel", event, cfg, utils, gt);
  mods_all.emplace_back(preselmod);
  ADDMOD(MapMod);
  if (analysis.unpackedGen)
    ADDMOD(DeepGenMod<UnpackedGenParticle>);
  else
    ADDMOD(DeepGenMod<GenParticle>);
  ADDMOD(TriggerMod);
  ADDMOD(SimpleLeptonMod);
  ADDMOD(ComplicatedLeptonMod);
  ADDMOD(SimplePhotonMod);
  ADDMOD(ComplicatedPhotonMod);
  ADDMOD(RecoilMod);
  ADDMOD(FatJetMod);
  ADDMOD(JetMod);
  ADDMOD(TauMod);

  postselmod = new ContainerMod("post-sel", event, cfg, utils, gt);
  mods_all.emplace_back(postselmod);
  ADDMOD(GenPMod);
  ADDMOD(JetFlavorMod); 
  ADDMOD(HbbMiscMod);
  ADDMOD(GenVHMod); 
  ADDMOD(KinFitMod);
  ADDMOD(InclusiveLeptonMod);
  ADDMOD(SoftActivityMod);
  ADDMOD(FatJetMatchingMod);
  ADDMOD(BTagSFMod);
  ADDMOD(BTagWeightMod);
  ADDMOD(TriggerEffMod);
  ADDMOD(GenStudyEWKMod);
  ADDMOD(QCDUncMod);
  ADDMOD(GenLepMod);
  ADDMOD(GenJetNuMod);
  ADDMOD(HFCountingMod);
  ADDMOD(KFactorMod);

  for (auto& mod : mods_all)
    mod->print();

  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Reading inputs");
  // Read inputs
  getInput();

  event.setStatus(*tIn, {"!*"});
  event.setAddress(*tIn, cfgmod.get_inputBranches());

  TH1D* hDTotalMCWeight = static_cast<TH1D*>(static_cast<TH1D*>(fIn->Get("hSumW"))->Clone("hDTotalMCWeight"));
  hDTotalMCWeight->SetDirectory(0);
  TH1D* hDNPUWeight = nullptr; {
    TH1D* hbase = static_cast<TH1D*>(fIn->Get("hNPVReco"));
    if (hbase == nullptr)
      hbase = static_cast<TH1D*>(fIn->Get("hNPVTrue"));
    if (hbase != nullptr) {
      hDNPUWeight = static_cast<TH1D*>(hbase->Clone("hDNPUWeight"));
      hDNPUWeight->SetDirectory(0);
    }
  }

  TTree* tW = static_cast<TTree*>(fIn->Get("weights"));
  if (tW && analysis.processType == kSignal) {
    if (tW->GetEntries()!=377 && tW->GetEntries()!=22) {
      logger.error("PandaAnalyzer::PandaAnalyzer",
          TString::Format("Reweighting failed because only found %u weights!",
                          unsigned(tW->GetEntries())));
      throw runtime_error("");
    }
    TString *id = new TString();
    tW->SetBranchAddress("id",&id);
    unsigned nW = tW->GetEntriesFast();
    for (unsigned iW=0; iW!=nW; ++iW) {
      tW->GetEntry(iW);
      wIDs->push_back(*id);
    }
  } else if (analysis.processType==kSignal) {
    logger.error("PandaAnalyzer::PandaAnalyzer","This is a signal file, but the weights are missing!");
    throw runtime_error("");
  }
  registry.publishConst("wIDs", wIDs);

  // Define outputs
  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Writing outputs");
  gt.is_monohiggs      = (analysis.monoh || analysis.hbb);
  gt.is_hbb            = analysis.hbb;
  gt.is_vh             = analysis.vqqhbb; 
  gt.is_vbf            = analysis.vbf;
  gt.is_fatjet         = (analysis.fatjet || analysis.deepGen);
  gt.is_leptonic       = analysis.complicatedLeptons;
  gt.is_photonic       = analysis.complicatedPhotons;
  gt.is_monotop        = !(analysis.monoh || analysis.hbb || analysis.vbf);
  gt.is_breg           = analysis.bjetRegTraining;
  gt.btagWeights       = analysis.btagWeights;
  gt.useCMVA           = analysis.useCMVA;
  for (auto& id : *wIDs)
    gt.signal_weights[id] = 1;

  makeOutput(); 

  fOut->WriteTObject(hDTotalMCWeight); delete hDTotalMCWeight; hDTotalMCWeight = nullptr;
  if (hDNPUWeight != nullptr) {
    fOut->WriteTObject(hDNPUWeight); delete hDNPUWeight; hDNPUWeight = nullptr;
  }

  event.rng.setSize(20);

<<<<<<< HEAD
  delete hDTotalMCWeight;
  if(hDNPUWeight) delete hDNPUWeight;
  
  delete bjetregReader;
  delete rochesterCorrection;
=======
  // read input data
  cfgmod.readData(analysis.datapath);
  for (auto& mod : mods_all)
    mod->readData(analysis.datapath);
>>>>>>> sid/master

  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Called constructor");
}


PandaAnalyzer::~PandaAnalyzer()
{
}

void PandaAnalyzer::AddGoodLumiRange(int run, int l0, int l1)
{
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) { // don't know about this run yet
    vector<LumiRange> newLumiList;
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
      logger.debug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i",run));
    return false;
  }

<<<<<<< HEAD
  jets = &event.chsAK4Jets;

  // these are bins of b-tagging eff in pT and eta, derived in 8024 TT MC
  // TODO: don't hardcode these 
  std::vector<double> vbtagpt {20.0,50.0,80.0,120.0,200.0,300.0,400.0,500.0,700.0,1000.0};
  std::vector<double> vbtageta {0.0,0.5,1.5,2.5};
  lfeff  = {{0.081,0.065,0.060,0.063,0.072,0.085,0.104,0.127,0.162},
            {0.116,0.097,0.092,0.099,0.112,0.138,0.166,0.185,0.222},
            {0.173,0.145,0.149,0.175,0.195,0.225,0.229,0.233,0.250}};
  ceff = {{0.377,0.389,0.391,0.390,0.391,0.375,0.372,0.392,0.435},
          {0.398,0.407,0.416,0.424,0.424,0.428,0.448,0.466,0.500},
          {0.375,0.389,0.400,0.425,0.437,0.459,0.481,0.534,0.488}};
  beff = {{0.791,0.815,0.825,0.835,0.821,0.799,0.784,0.767,0.760},
          {0.794,0.816,0.829,0.836,0.823,0.804,0.798,0.792,0.789},
          {0.739,0.767,0.780,0.789,0.776,0.771,0.779,0.787,0.806}};
  btagpt = Binner(vbtagpt);
  btageta = Binner(vbtageta);

  std::vector<unsigned int> metTriggers;
  std::vector<unsigned int> eleTriggers;
  std::vector<unsigned int> phoTriggers;
  std::vector<unsigned int> muTriggers;
  std::vector<unsigned int> jetTriggers;
  std::vector<unsigned int> muFakeTriggers;
  std::vector<unsigned int> eleFakeTriggers;

  if (isData || analysis->applyMCTriggers) {
    if (DEBUG) PDebug("PandaAnalyzer::Run","Loading the trigger paths");
    std::vector<TString> paths;
    paths = {
          "HLT_PFMET170_NoiseCleaned",
          "HLT_PFMET170_HBHECleaned",
          "HLT_PFMET170_JetIdCleaned",
          "HLT_PFMET170_NotCleaned",
          "HLT_PFMET170_HBHE_BeamHaloCleaned",
          "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight",
          "HLT_PFMETNoMu110_NoiseCleaned_PFMHTNoMu110_IDTight",
          "HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight",
          "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight",
          "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight",
          "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
          "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"
    };
    triggerHandlers[kMETTrig].addTriggers(paths);

    if (analysis->complicatedLeptons)
      paths = {
# 2016
          "HLT_Ele25_eta2p1_WPTight_Gsf",
          "HLT_Ele27_eta2p1_WPLoose_Gsf",
          "HLT_Ele27_WPTight_Gsf",
          "HLT_Ele30_WPTight_Gsf",
          "HLT_Ele35_WPLoose_Gsf",
          "HLT_Ele27_WP85_Gsf",
          "HLT_Ele27_WPLoose_Gsf",
          "HLT_Ele105_CaloIdVT_GsfTrkIdT",
          "HLT_Ele115_CaloIdVT_GsfTrkIdT",
          "HLT_Ele27_eta2p1_WPTight_Gsf",
          "HLT_Ele32_eta2p1_WPTight_Gsf",
          "HLT_ECALHT800",
# 2017
          "HLT_Ele115_CaloIdVT_GsfTrkIdT",
	  "HLT_Ele27_WPTight_Gsf",
	  "HLT_Ele32_WPTight_Gsf",
	  "HLT_Ele35_WPTight_Gsf",
	  "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
	  "HLT_Photon200"
      };
    else
      paths = {
            "HLT_Ele27_WP85_Gsf",
            "HLT_Ele27_WPLoose_Gsf",
            "HLT_Ele105_CaloIdVT_GsfTrkIdT",
            "HLT_Ele27_WPTight_Gsf",
            "HLT_Ele30_WPTight_Gsf",
            "HLT_Ele27_eta2p1_WPTight_Gsf",
            "HLT_Ele32_eta2p1_WPTight_Gsf",
            "HLT_Ele35_WPLoose_Gsf",
            "HLT_ECALHT800",
	  
	  "HLT_Ele32_WPTight_Gsf"
      };
    triggerHandlers[kSingleEleTrig].addTriggers(paths);
    
    paths = {
# 2016
          "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
          "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
          "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
          "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
# 2017	  
	  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
	  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
	  "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8",
	  "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8"
    };
    triggerHandlers[kDoubleMuTrig].addTriggers(paths);

    paths = {
# 2016
          "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
          "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf",
# 2017	  	  
	  "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
	  "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
	  "HLT_DiEle27_WPTightCaloOnly_L1DoubleEG",
	  "HLT_DoubleEle33_CaloIdL_MW",
	  "HLT_DoubleEle25_CaloIdL_MW",
	  "HLT_DoublePhoton70"
    };
    triggerHandlers[kDoubleEleTrig].addTriggers(paths);
    
    paths = {
# 2016
          "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
          "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
          "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
          "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
          "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
          "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
          "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
          "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
# 2017
          "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
          "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
          "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
          "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
          "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
          "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
          "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
          "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"
    };
    triggerHandlers[kEMuTrig].addTriggers(paths);

    paths = {
          "HLT_Photon175",
          "HLT_Photon165_HE10",
          "HLT_Photon36_R9Id90_HE10_IsoM",
          "HLT_Photon50_R9Id90_HE10_IsoM",
          "HLT_Photon75_R9Id90_HE10_IsoM",
          "HLT_Photon90_R9Id90_HE10_IsoM",
          "HLT_Photon120_R9Id90_HE10_IsoM",
          "HLT_Photon165_R9Id90_HE10_IsoM",
          "HLT_Photon300_NoHE",
          "HLT_ECALHT800",
	  
	  "HLT_Photon200"
    };
    triggerHandlers[kSinglePhoTrig].addTriggers(paths);

    paths = {
          "HLT_PFHT650",
          "HLT_PFHT900",
          "HLT_PFJet500",
          "HLT_PFJet450",
          "HLT_PFJet320",
    };
    triggerHandlers[kJetHTTrig].addTriggers(paths);

    if (analysis->complicatedLeptons)
      paths = {
# 2016
          "HLT_IsoMu24",
          "HLT_IsoTkMu24",
          "HLT_IsoMu22",
          "HLT_IsoTkMu22",
          "HLT_Mu45_eta2p1",
          "HLT_Mu50",
# 2017
          "HLT_IsoMu24",	  
	  "HLT_IsoMu27",
	  "HLT_IsoMu30"
          "HLT_Mu50",
	  "HLT_TkMu100"
      };
    else
      paths = {
            "HLT_IsoMu20",
            "HLT_IsoMu22",
            "HLT_IsoMu24",

	    "HLT_IsoMu27"
      };
    triggerHandlers[kSingleMuTrig].addTriggers(paths);
    
    paths = {
# 2016
          "HLT_Mu8_TrkIsoVVL",
          "HLT_Mu17_TrkIsoVVL",
# 2017
          "HLT_Mu8_TrkIsoVVL",
          "HLT_Mu17_TrkIsoVVL",	  
	  "HLT_Mu19_TrkIsoVVL"
    };
    triggerHandlers[kMuFakeTrig].addTriggers(paths);

    paths = {
# 2016
          "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",
          "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",
          "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30",
          "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30"
# 2017
          "HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",
          "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",
          "HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30",
          "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30"
    };
    triggerHandlers[kEleFakeTrig].addTriggers(paths);

    RegisterTriggers();
=======
  // found the run, now look for a lumi range
  for (auto &range : run_->second) {
    if (range.Contains(lumi)) {
      if (DEBUG)
        logger.debug("PandaAnalyzer::PassGoodLumis",TString::Format("Accepting run=%i, lumi=%i",run,lumi));
      return true;
    }
>>>>>>> sid/master
  }

  // matched no lumi range
  if (DEBUG)
    logger.debug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i, lumi=%i",run,lumi));
  return false;
}


bool PandaAnalyzer::PassPresel(Selection::Stage stage)
{
  if (selections.size() == 0)
    return true;

  bool pass = false;
  for (auto& s : selections) {
    if (s->anded())
      continue;
    if (DEBUG>1)
      logger.debug("PandaAnalyzer::PassPresel",s->get_name());
    if (s->accept(stage)) {
      pass = true;
      break;
    }
  }

  for (auto& s : selections) {
    if (s->anded()) {
      if (DEBUG>1)
        logger.debug("PandaAnalyzer::PassPresel",s->get_name());
      pass = pass && s->accept(stage);
    }
  }

  return pass;
}

<<<<<<< HEAD
    tr->TriggerEvent(TString::Format("GetEntry %u",iE));
    if (DEBUG > 5) {
      PDebug("PandaAnalyzer::Run::Dump","");
      event.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.photons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.muons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.electrons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.chsAK4Jets.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      event.pfMet.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaAnalyzer::Run::Dump","");
      //event.metMuOnlyFix.print(std::cout, 2);
      std::cout << std::endl;
    }
    
    gt->filter_maxRecoil = event.recoil.max;
=======
>>>>>>> sid/master


void PandaAnalyzer::Reset()
{
  for (auto& mod : mods_all)
    mod->reset();

  Analyzer::Reset();
}



void PandaAnalyzer::Terminate()
{
  for (auto& mod : mods_all)
    mod->terminate();

  Analyzer::Terminate();
}


// run
void PandaAnalyzer::Run()
{

  // INITIALIZE --------------------------------------------------------------------------
  unsigned nZero, nEvents, iE=0;
  setupRun(nZero, nEvents); 

  for (auto& mod : mods_all)
    mod->initialize(registry);

  ProgressReporter pr("PandaAnalyzer::Run",&iE,&nEvents,100);
  TimeReporter& tr = cfgmod.cfg.tr;
  tr.TriggerEvent("configuration"); 

  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    pr.Report();

    Reset();
    event.getEntry(*tIn,iE);
    tr.TriggerEvent(TString::Format("GetEntry %u",iE));

    gblmod->execute();

    if (analysis.isData && !PassGoodLumis(gt.runNumber,gt.lumiNumber))
        continue;

    preselmod->execute();

    if (!PassPresel(Selection::sReco))
      continue;

    postselmod->execute();

    if (!PassPresel(Selection::sGen)) // only check gen presel here
      continue;

//     if (analysis.deep || analysis.hbb)
//       FatJetPartons();
//     if (analysis.deep) {
//       FillPFTree();
//       tAux->Fill();
//       if (tAux->GetEntriesFast() == 2500)
//         IncrementAuxFile();
//       tr->TriggerEvent("aux fill");
//     }

    gt.Fill();

    tr.TriggerEvent("fill");

  }

  tr.Summary();
  for (auto& s : selections)
    s->report();

  if (DEBUG) { logger.debug("PandaAnalyzer::Run","Done with entry loop"); }

} // Run()
