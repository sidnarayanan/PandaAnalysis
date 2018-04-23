#include "../interface/PandaAnalyzer.h"
#include "TMath.h"
#include <algorithm>
#include <vector>


// only contain the main methods

using namespace panda;
using namespace std;

PandaAnalyzer::PandaAnalyzer(int debug_/*=0*/) 
{
  DEBUG = debug_;

  if (DEBUG) PDebug("PandaAnalyzer::PandaAnalyzer","Calling constructor");
  gt = new GeneralTree();
  if (DEBUG) PDebug("PandaAnalyzer::PandaAnalyzer","Built GeneralTree");
  ibetas = gt->get_ibetas();
  Ns = gt->get_Ns();
  orders = gt->get_orders();
  if (DEBUG) PDebug("PandaAnalyzer::PandaAnalyzer","Called constructor");
}


PandaAnalyzer::~PandaAnalyzer() 
{
  if (DEBUG) PDebug("PandaAnalyzer::~PandaAnalyzer","Calling destructor");
}


void PandaAnalyzer::ResetBranches() 
{
  genObjects.clear();
  matchPhos.clear();
  matchEles.clear();
  matchLeps.clear();
  looseLeps.clear();
  tightLeps.clear();
  loosePhos.clear();
  genJetsNu.clear();
  validGenP.clear();
  for (int i = 0; i != jes2i(shiftjes::N); ++i) 
    jesShifts[i].clear();
  fj1 = 0;
  gt->Reset();
  if (DEBUG) PDebug("PandaAnalyzer::ResetBranches","Reset");
}



void PandaAnalyzer::Terminate() 
{
  fOut->WriteTObject(tOut);
  fOut->Close();
  fOut = 0; tOut = 0;

  if (analysis->deep)
    IncrementAuxFile(true);
  if (analysis->deepGen)
    IncrementGenAuxFile(true);

  for (unsigned i = 0; i != cN; ++i) {
    delete h1Corrs[i];
    h1Corrs[i] = 0;
  }
  for (unsigned i = 0; i != cN; ++i) {
    delete h2Corrs[i];
    h2Corrs[i] = 0;
  }
  for (auto *f : fCorrs)
    if (f)
      f->Close();

  delete btagCalib;
  delete sj_btagCalib;
  for (auto *reader : btagReaders)
    delete reader;

  for (auto& iter : ak8UncReader)
    delete iter.second;

  delete ak8JERReader;

  for (auto& iter : ak4UncReader)
    delete iter.second;

  for (auto& iter : ak4ScaleReader) {
    delete iter.second;
  }

  delete ak4JERReader;

  delete activeArea;
  delete areaDef;
  delete jetDef;
  delete jetDefKt;
  delete jetDefGen;
  delete softDrop;

  delete hDTotalMCWeight;
  
  delete bjetregReader;
  delete rochesterCorrection;

  delete ecfcalc;
  delete grid;

  if (DEBUG) PDebug("PandaAnalyzer::Terminate","Finished with output");
}

void shiftMET(const panda::RecoMet& met, TLorentzVector& v, shiftjes shift) 
{
  float pt;
  float phi;
  switch (shift) {
    case shiftjes::kNominal:
      pt = met.pt;
      phi = met.phi;
      break;
    case shiftjes::kJESUp:
      pt = met.ptCorrUp;
      phi = met.phiCorrUp;
      break;
    case shiftjes::kJESDown:
      pt = met.ptCorrDown;
      phi = met.phiCorrDown;
      break;
    default:
      PError("shiftMET", "Unknown JES type!");
      exit(1);
  }

  v.SetPtEtaPhiM(pt, 0, phi, 0);

}


// run
void PandaAnalyzer::Run() 
{

  fOut->cd(); // to be absolutely sure

  // INITIALIZE --------------------------------------------------------------------------

  unsigned int nEvents = tIn->GetEntries();
  unsigned int nZero = 0;
  if (lastEvent>=0 && lastEvent<(int)nEvents)
    nEvents = lastEvent;
  if (firstEvent>=0)
    nZero = firstEvent;

  if (!fOut || !tIn) {
    PError("PandaAnalyzer::Run","NOT SETUP CORRECTLY");
    exit(1);
  }

  // get bounds
  genBosonPtMin=150, genBosonPtMax=1000;
  if (!isData && h1Corrs[cZNLO]) {
    genBosonPtMin = h1Corrs[cZNLO]->GetHist()->GetBinCenter(1);
    genBosonPtMax = h1Corrs[cZNLO]->GetHist()->GetBinCenter(h1Corrs[cZNLO]->GetHist()->GetNbinsX());
  }

  if (analysis->ak8) {
    if (analysis->puppiJets)
      fatjets = &event.puppiAK8Jets;
    else
      fatjets = &event.chsAK8Jets;
  } else if (analysis->fatjet) {
    if (analysis->puppiJets)
      fatjets = &event.puppiCA15Jets;
    else
      fatjets = &event.chsCA15Jets;
  }

  ak4jets = &event.chsAK4Jets;

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
          "HLT_ECALHT800"
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
            "HLT_ECALHT800"
      };
    triggerHandlers[kSingleEleTrig].addTriggers(paths);
    
    paths = {
          "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
          "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
          "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
          "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"
    };
    triggerHandlers[kDoubleMuTrig].addTriggers(paths);

    paths = {
          "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
          "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf"
    };
    triggerHandlers[kDoubleEleTrig].addTriggers(paths);
    
    paths = {
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
          "HLT_ECALHT800"
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
          "HLT_IsoMu24",
          "HLT_IsoTkMu24",
          "HLT_IsoMu22",
          "HLT_IsoTkMu22",
          "HLT_Mu45_eta2p1",
          "HLT_Mu50"
      };
    else
      paths = {
            "HLT_IsoMu20",
            "HLT_IsoMu22",
            "HLT_IsoMu24",
      };
    triggerHandlers[kSingleMuTrig].addTriggers(paths);
    
    paths = {
          "HLT_Mu8_TrkIsoVV",
          "HLT_Mu17_TrkIsoVV"
    };
    triggerHandlers[kMuFakeTrig].addTriggers(paths);

    paths = {
          "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",
          "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30",
          "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30"
    };
    triggerHandlers[kEleFakeTrig].addTriggers(paths);

    RegisterTriggers();
  }

  if (analysis->ak8)
    FATJETMATCHDR2 = 0.64;

  fOut->cd(); // to be absolutely sure

  // set up reporters
  unsigned int iE=0;
  ProgressReporter pr("PandaAnalyzer::Run",&iE,&nEvents,10);
  tr = new TimeReporter("PandaAnalyzer::Run",DEBUG+1);


  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    tr->Start();
    pr.Report();
    ResetBranches();
    event.getEntry(*tIn,iE);


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
      event.metMuOnlyFix.print(std::cout, 2);
      std::cout << std::endl;
    }
    
    gt->filter_maxRecoil = event.recoil.max;

    if (!PassPresel(Selection::sRecoil))
      continue;

    // event info
    gt->mcWeight = event.weight;
    gt->runNumber = event.runNumber;
    gt->lumiNumber = event.lumiNumber;
    gt->eventNumber = event.eventNumber;
    gt->isData = isData ?  1 : 0; 
    gt->npv = event.npv;
    gt->pu = event.npvTrue;
    gt->metFilter = (event.metFilters.pass()) ? 1 : 0;
    gt->metFilter = (gt->metFilter==1 && !event.metFilters.badPFMuons) ? 1 : 0;
    gt->metFilter = (gt->metFilter==1 && !event.metFilters.badChargedHadrons) ? 1 : 0;

    if (isData) {
      // check the json
      if (!PassGoodLumis(gt->runNumber,gt->lumiNumber))
        continue;

    } else { // !isData
      gt->sf_npv = GetCorr(cNPV,gt->npv);
      gt->sf_pu = GetCorr(cPU,gt->pu);
    }

    // save triggers
    if (isData || analysis->applyMCTriggers) {
      for (unsigned iT = 0; iT != kNTrig; ++iT) {
        auto &th = triggerHandlers.at(iT);
        for (auto iP : th.indices) {
          if (event.triggerFired(iP)) {
              gt->trigger |= (1 << iT);
              break;
          }
        }
      }
    }

    if (analysis->rerunJES)
      SetupJES();

    tr->TriggerEvent("initialize");

    // met
    gt->pfmetRaw = event.rawMet.pt;
    gt->calomet = event.caloMet.pt;
    gt->sumETRaw = event.pfMet.sumETRaw;
    gt->trkmet = event.trkMet.pt;
    gt->trkmetphi = event.trkMet.phi;

    JESLOOP {
      auto& jets = jesShifts[shift];
      // PF 
      shiftMET(event.pfMet, jets.vpfMET, i2jes(shift));
      gt->pfmet[shift] = jets.vpfMET.Pt();
      gt->pfmetphi[shift] = jets.vpfMET.Phi();

      // Puppi
      shiftMET(event.puppiMet, jets.vpuppiMET, i2jes(shift));
      gt->puppimet[shift] = jets.vpuppiMET.Pt();
      gt->puppimetphi[shift] = jets.vpuppiMET.Phi();

      jets.vpfMETNoMu.SetMagPhi(gt->pfmet[shift], gt->pfmetphi[shift]);
    }


    tr->TriggerEvent("met");

    if (!isData) {
      // take care of bug in panda version <= 009
      // replace duplicate gen particles with preferred version.
      // template type could be inferred but let's be explicit
      // so we can read it.
      if (event.genParticles.size() > 0) {
        RemoveGenDups<GenParticle>(event.genParticles);
      } else {
        RemoveGenDups<UnpackedGenParticle>(event.genParticlesU);
      }

      // do this up here before the preselection
      if (analysis->deepGen) {
        if (event.genParticles.size() > 0) 
          FillGenTree<GenParticle>();
        else
          FillGenTree<UnpackedGenParticle>();
        if (gt->genFatJetPt > 400 || analysis->deepExC) 
          tAux->Fill();
        if (tAux->GetEntriesFast() == 2500)
          IncrementGenAuxFile();
        tr->TriggerEvent("fill gen aux");
      }
    }


    if (!analysis->genOnly) {
      // electrons and muons
      if (analysis->complicatedLeptons) {
        ComplicatedLeptons();
      } else {
        SimpleLeptons();
      }
      
      // photons
      if (analysis->complicatedPhotons) {
        ComplicatedPhotons();
      } else {
        SimplePhotons();
      }

      // recoil!
      if (analysis->recoil)
        Recoil();

      // fatjets
      if (analysis->fatjet) {
        FatjetBasics();
        if (analysis->recluster)
          FatjetRecluster();
        tr->TriggerEvent("fatjet");
      }

      // interesting jets
      JetBasics();

      Taus();

      if (!PassPresel(Selection::sReco)) // only check reco presel here
        continue;

      if (analysis->hbb) {
        JetHbbSoftActivity();
        GetMETSignificance();
      }
    }

    if (!isData) {
      if (!analysis->genOnly) {
        if (analysis->fatjet)
          FatjetMatching();

        if (analysis->btagSFs)
          JetBtagSFs();
        if (analysis->btagWeights)
          JetCMVAWeights();
        
        TriggerEffs();

        if (analysis->complicatedLeptons ||
            analysis->complicatedPhotons)
          GenStudyEWK();
        else
          LeptonSFs();

        PhotonSFs();
      }

      QCDUncs();
      SignalReweights();

      if (analysis->vbf)
        SaveGenLeptons();

      SignalInfo();

      if (analysis->reclusterGen && analysis->hbb) {
        GenJetsNu();
        MatchGenJets(genJetsNu);
      }

      if (analysis->hfCounting)
        HeavyFlavorCounting();

      TopPTReweight();
      VJetsReweight();
    }

    
    if (!PassPresel(Selection::sGen)) // only check gen presel here
      continue;

    if (analysis->deep) {
      FatjetPartons();
      FillPFTree();
      tAux->Fill();
      if (tAux->GetEntriesFast() == 2500)
        IncrementAuxFile();
      tr->TriggerEvent("aux fill");
    }

    gt->Fill();

    tr->TriggerEvent("fill");

  } // entry loop

  tr->Summary();
  for (auto* s : selections) 
    s->report(); 

  if (DEBUG) { PDebug("PandaAnalyzer::Run","Done with entry loop"); }

} // Run()

