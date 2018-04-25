#include "../interface/CommonMods.h"

using namespace pa;
using namespace std;
using namespace panda;

void TriggerMod::do_init(Registry& registry) 
{
  vector<TString> paths; 
  if (analysis.isData || analysis.applyMCTriggers) {
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

    for (auto &th : triggerHandlers) {
      unsigned N = th.paths.size();
      for (unsigned i = 0; i != N; i++) {
        unsigned panda_idx = event.registerTrigger(th.paths.at(i));
        th.indices[i] = panda_idx;
      }
    }
  }
}

void TriggerMod::do_execute()
{
  if (analysis.isData || analysis.applyMCTriggers) {
    for (unsigned iT = 0; iT != kNTrig; ++iT) {
      auto &th = triggerHandlers.at(iT);
      for (auto iP : th.indices) {
        if (event.triggerFired(iP)) {
            gt.trigger |= (1 << iT);
            break;
        }
      }
    }
  }
}

void GlobalMod::do_init(Registry& registry)
{
  registry.publish("jesShifts", &jesShifts);
}

void GlobalMod::do_execute()
{
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
  
  gt.filter_maxRecoil = event.recoil.max;

  // event info
  gt.mcWeight = event.weight;
  gt.runNumber = event.runNumber;
  gt.lumiNumber = event.lumiNumber;
  gt.eventNumber = event.eventNumber;
  gt.isData = isData ?  1 : 0; 
  gt.npv = event.npv;
  gt.pu = event.npvTrue;
  gt.metFilter = (event.metFilters.pass()) ? 1 : 0;
  gt.metFilter = (gt.metFilter==1 && !event.metFilters.badPFMuons) ? 1 : 0;
  gt.metFilter = (gt.metFilter==1 && !event.metFilters.badChargedHadrons) ? 1 : 0;

  if (!analysis.isData) {
    gt.sf_npv = utils.getCorr(cNPV, gt.npv);
    gt.sf_pu = utils.getCorr(cPU, gt.pu);
  }

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
}
