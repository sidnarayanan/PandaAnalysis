#include "../interface/CommonMods.h"

using namespace pa;
using namespace std;
using namespace panda;

void RecoilMod::do_execute()
{
  TLorentzVector vObj1, vObj2;
  gt.whichRecoil = 0; // -1=photon, 0=MET, 1,2=nLep
  if (gt.nLooseLep>0) {
    panda::Lepton *lep1 = looseLeps->at(0);
    vObj1.SetPtEtaPhiM(lep1->pt(),lep1->eta(),lep1->phi(),lep1->m());

    // one lep => W
    JESLOOP {
      auto& jets = jesShifts[shift]; 

      jets.vpuppiUW = jets.vpuppiMET + vObj1;
      gt.puppiUWmag[shift] = jets.vpuppiUW.Pt(); 
      gt.puppiUWphi[shift] = jets.vpuppiUW.Phi(); 

      jets.vpfUW = jets.vpfMET + vObj1;
      gt.pfUWmag[shift] = jets.vpfUW.Pt(); 
      gt.pfUWphi[shift] = jets.vpfUW.Phi(); 
    }

    if (gt.nLooseLep>1 && lepPdgId[0]+lepPdgId[1]==0) {
      // two OS lep => Z
      panda::Lepton *lep2 = looseLeps->at(1);
      vObj2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());

      JESLOOP {
        auto& jets = jesShifts[shift]; 

        jets.vpuppiUZ = jets.vpuppiUW + vObj2;
        gt.puppiUZmag[shift] = jets.vpuppiUZ.Pt(); 
        gt.puppiUZphi[shift] = jets.vpuppiUZ.Phi(); 

        jets.vpfUZ = jets.vpfUW + vObj2;
        gt.pfUZmag[shift] = jets.vpfUZ.Pt(); 
        gt.pfUZphi[shift] = jets.vpfUZ.Phi(); 
      }

      gt.whichRecoil = 2;
    } else {
      gt.whichRecoil = 1;
    }
  }
  if (gt.nLoosePhoton>0) {
    panda::Photon *pho = loosePhos->(0);
    vObj1.SetPtEtaPhiM(pho->pt(),pho->eta(),pho->phi(),0.);

    JESLOOP {
      auto& jets = jesShifts[shift]; 

      jets.vpuppiUA = jets.vpuppiMET + vObj1;
      gt.puppiUAmag[shift] = jets.vpuppiUA.Pt(); 
      gt.puppiUAphi[shift] = jets.vpuppiUA.Phi(); 

      jets.vpfUA = jets.vpfMET + vObj1;
      gt.pfUAmag[shift] = jets.vpfUA.Pt(); 
      gt.pfUAphi[shift] = jets.vpfUA.Phi(); 
    }

    if (gt.nLooseLep==0) {
      gt.whichRecoil = -1;
    }
  }
  if (gt.nLooseLep==0 && gt.nLoosePhoton==0) {
    gt.whichRecoil = 0;
  }
}

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

void GlobalMod::do_execute()
{
  if (cfg.DEBUG > 5) {
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

  gt.pfmetRaw = event.rawMet.pt;
  gt.calomet = event.caloMet.pt;
  gt.sumETRaw = event.pfMet.sumETRaw;
  gt.trkmet = event.trkMet.pt;
  gt.trkmetphi = event.trkMet.phi;
  
  JESLOOP {
    auto& jets = jesShifts[shift];
    // PF 
    shiftMET(event.pfMet, jets.vpfMET, i2jes(shift));
    gt.pfmet[shift] = jets.vpfMET.Pt();
    gt.pfmetphi[shift] = jets.vpfMET.Phi();

    // Puppi
    shiftMET(event.puppiMet, jets.vpuppiMET, i2jes(shift));
    gt.puppimet[shift] = jets.vpuppiMET.Pt();
    gt.puppimetphi[shift] = jets.vpuppiMET.Phi();

    jets.vpfMETNoMu.SetMagPhi(gt.pfmet[shift], gt.pfmetphi[shift]);
  }
}

void GenPMod::do_execute()
{
  if (analysis.isData)
    return;
  if (event.genParticles.size() > 0) {
    merge_particles(event.genParticles);
  } else {
    merge_particles(event.genParticlesU);
  }
}

