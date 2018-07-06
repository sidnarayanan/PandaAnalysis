#include "../interface/CommonMods.h"

using namespace pa;
using namespace std;
using namespace panda;

void shiftMET(const RecoMet& met, TLorentzVector& v, shiftjes shift) 
{
  float pt;
  float phi;
  switch (shift) {
    case shiftjes::kNominal:
      pt = met.pt;
      phi = met.phi;
      break;
    case shiftjes::kJESTotalUp:
      pt = met.ptCorrUp;
      phi = met.phiCorrUp;
      break;
    case shiftjes::kJESTotalDown:
      pt = met.ptCorrDown;
      phi = met.phiCorrDown;
      break;
    default:
      logger.error("shiftMET", "Unknown JES type!");
      exit(1);
  }

  v.SetPtEtaPhiM(pt, 0, phi, 0);
}

void TriggerMod::do_init(Registry& registry) 
{
  vector<TString> paths; 
  if (analysis.isData || analysis.mcTriggers) {
    // MET
    if (analysis.year == 2016) {
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
    } else if (analysis.year == 2017) { 
        paths = {
          "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
          "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
          "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",
          "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",
        };
    }
    triggerHandlers[kMETTrig].addTriggers(paths);

    // SingleEle    
    if (analysis.complicatedLeptons) {
      if (analysis.year == 2016) {
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
      } else if (analysis.year == 2017) {
        paths = {
          "HLT_Ele115_CaloIdVT_GsfTrkIdT",
          "HLT_Ele27_WPTight_Gsf",
          "HLT_Ele32_WPTight_Gsf",
          "HLT_Ele35_WPTight_Gsf",
          "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
          "HLT_Photon200"
        };
      }
    } else {
      if (analysis.year == 2016) {
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
      } else if (analysis.year == 2017) {
        paths = {
          "HLT_Ele35_WPTight_Gsf",
          "HLT_Ele38_WPTight_Gsf",
          "HLT_Ele40_WPTight_Gsf"
        };
      }
    }
    triggerHandlers[kSingleEleTrig].addTriggers(paths);
    
    // single pho
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
          "HLT_CaloJet500_NoJetID"
    };
    triggerHandlers[kSinglePhoTrig].addTriggers(paths);

    // Single muon
    if (analysis.complicatedLeptons || analysis.recalcECF) { // either comp lepton or tnp
      if (analysis.year == 2016) {
        paths = {
          "HLT_IsoMu24",
          "HLT_IsoTkMu24",
          "HLT_IsoMu22",
          "HLT_IsoTkMu22",
          "HLT_Mu45_eta2p1",
          "HLT_Mu50"
        };
      } else if (analysis.year == 2017) {
        paths = {
          "HLT_IsoMu27"
          "HLT_IsoMu24",      
          "HLT_IsoMu27",
          "HLT_IsoMu30"
          "HLT_Mu50"
        };
      }
    } else {
      if (analysis.year == 2016) {
        paths = {
          "HLT_IsoMu20",
          "HLT_IsoMu22",
          "HLT_IsoMu24",
        };
      } else if (analysis.year == 2017) {
        paths = {
          "HLT_IsoMu24",
          "HLT_IsoMu27"
        };
      }
    }
    triggerHandlers[kSingleMuTrig].addTriggers(paths);

    // double muon
    if (analysis.year==2016) { 
      paths = {
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"
      };
    } else if (analysis.year==2017) {
      paths = {
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8",
        "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8"
      };
    }
    triggerHandlers[kDoubleMuTrig].addTriggers(paths);

    // double ele
    if (analysis.year==2016) 
      paths = {
            "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf"
      };
    else if (analysis.year==2017)
      paths = {
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
        "HLT_DiEle27_WPTightCaloOnly_L1DoubleEG",
        "HLT_DoubleEle33_CaloIdL_MW",
        "HLT_DoubleEle25_CaloIdL_MW",
        "HLT_DoublePhoton70"
      };
    triggerHandlers[kDoubleEleTrig].addTriggers(paths);
    
    // emu
    if (analysis.year==2016) {
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
    } else if (analysis.year==2017) {
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
    }
    triggerHandlers[kEMuTrig].addTriggers(paths);

    // JetHT
    paths = {
          "HLT_PFHT650",
          "HLT_PFHT900",
          "HLT_PFJet500",
          "HLT_PFJet450",
          "HLT_PFJet320",
    };
    triggerHandlers[kJetHTTrig].addTriggers(paths);

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
    if (analysis.hbb && analysis.year == 2017) {
      event.registerTriggerObjects("hltEle32L1DoubleEGWPTightGsfTrackIsoFilter");
      event.registerTriggerObjects("hltEGL1SingleEGOrFilter");
    }
  }
}

void TriggerMod::checkEle32()
{
  auto& filter1Objects = event.triggerObjects.filterObjects("hltEle32L1DoubleEGWPTightGsfTrackIsoFilter");
  auto& filter2Objects = event.triggerObjects.filterObjects("hltEGL1SingleEGOrFilter");
  if (filter1Objects.size()==0 || filter2Objects.size()==0) 
    return; 
  TLorentzVector filter1ObjectP4 = filter1Objects[0]->p4();
  TLorentzVector filter2ObjectP4 = filter2Objects[0]->p4();
  bool matchedToTriggerObject = false;
  for (auto& ele : event.electrons) {
    if (!ele.tight && !ele.mvaWP80) 
      continue;
    if (filter1ObjectP4.DeltaR(ele.p4())<0.1 && filter2ObjectP4.DeltaR(ele.p4())<0.1) {
      matchedToTriggerObject = true;
      break;
    }
  }
  if (matchedToTriggerObject) { 
    gt.trigger |= (1 << kSingleEleTrig);
  }
}

void TriggerMod::do_execute()
{
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

void GlobalMod::do_execute()
{
  if (cfg.DEBUG > 5) {
    logger.debug("PandaAnalyzer::Run::Dump","");
    event.print(std::cout, 2);
    std::cout << std::endl;
    logger.debug("PandaAnalyzer::Run::Dump","");
    event.photons.print(std::cout, 2);
    std::cout << std::endl;
    logger.debug("PandaAnalyzer::Run::Dump","");
    event.muons.print(std::cout, 2);
    std::cout << std::endl;
    logger.debug("PandaAnalyzer::Run::Dump","");
    event.electrons.print(std::cout, 2);
    std::cout << std::endl;
    logger.debug("PandaAnalyzer::Run::Dump","");
    event.chsAK4Jets.print(std::cout, 2);
    std::cout << std::endl;
    logger.debug("PandaAnalyzer::Run::Dump","");
    event.pfMet.print(std::cout, 2);
    std::cout << std::endl;
  }
  
  gt.filter_maxRecoil = event.recoil.max;

  // event info
  gt.mcWeight = event.weight;
  gt.runNumber = event.runNumber;
  gt.lumiNumber = event.lumiNumber;
  gt.eventNumber = event.eventNumber;
  gt.isData = analysis.isData ?  1 : 0; 
  gt.npv = event.npv;
  gt.rho = event.rho;
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
  gt.pfmetsig = event.pfMet.significance;
  gt.puppimetsig = event.puppiMet.significance; 
  
  METLOOP {
    auto& jets = (*jesShifts)[shift];
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

template <typename TREE>
void BaseGenPMod<TREE>::do_execute()
{
  if (this->event.genParticles.size() > 0) {
    merge_particles(this->event.genParticles);
  } else {
    merge_particles(this->event.genParticlesU);
  }
}

template class BaseGenPMod<GeneralTree>;
template class BaseGenPMod<HeavyResTree>;

void RecoilMod::do_execute()
{
  TLorentzVector vObj1, vObj2;
  gt.whichRecoil = 0; // -1=photon, 0=MET, 1,2=nLep
  if (gt.nLooseLep>0) {
    Lepton *lep1 = looseLeps->at(0);
    vObj1.SetPtEtaPhiM(lep1->pt(),lep1->eta(),lep1->phi(),lep1->m());

    // one lep => W
    METLOOP {
      auto& jets = (*jesShifts)[shift]; 

      jets.vpuppiUW = jets.vpuppiMET + vObj1;
      gt.puppiUWmag[shift] = jets.vpuppiUW.Pt(); 
      gt.puppiUWphi[shift] = jets.vpuppiUW.Phi(); 

      jets.vpfUW = jets.vpfMET + vObj1;
      gt.pfUWmag[shift] = jets.vpfUW.Pt(); 
      gt.pfUWphi[shift] = jets.vpfUW.Phi(); 
    }

    if (gt.nLooseLep>1 && (*lepPdgId)[0]+(*lepPdgId)[1]==0) {
      // two OS lep => Z
      Lepton *lep2 = looseLeps->at(1);
      vObj2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());

      METLOOP {
        auto& jets = (*jesShifts)[shift]; 

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
    Photon *pho = (*loosePhos)[0];
    vObj1.SetPtEtaPhiM(pho->pt(),pho->eta(),pho->phi(),0.);

    METLOOP {
      auto& jets = (*jesShifts)[shift]; 

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

void TriggerEffMod::do_execute() 
{
  // trigger efficiencies
  gt.sf_metTrig = utils.getCorr(cTrigMET,gt.pfmetnomu[jes2i(shiftjes::kNominal)]);
  gt.sf_metTrigZmm = utils.getCorr(cTrigMETZmm,gt.pfmetnomu[jes2i(shiftjes::kNominal)]);

  auto* lep0 = looseLeps->size()>0 ? (*looseLeps)[0] : nullptr;
  auto* lep1 = looseLeps->size()>1 ? (*looseLeps)[1] : nullptr;

  if (gt.nLooseElectron>1 && analysis.complicatedLeptons) {
    gt.sf_eleTrig = utils.getCorr(cTrigDoubleEleLeg1, gt.electronEta[0], gt.electronPt[0]) *
                    utils.getCorr(cTrigDoubleEleLeg2, gt.electronEta[1], gt.electronPt[1]);
  } else if (gt.nLooseElectron>0) {
    Electron *ele1=nullptr, *ele2=nullptr;
    if (gt.nLooseLep>0) ele1 = dynamic_cast<Electron*>(lep0);
    if (gt.nLooseLep>1) ele2 = dynamic_cast<Electron*>(lep1);
    float eff1=0, eff2=0;
    if (ele1!=nullptr && ele1->tight) {
      eff1 = utils.getCorr(cTrigEle, ele1->eta(), ele1->pt());
      if (ele2!=nullptr && ele2->tight)
        eff2 = utils.getCorr(cTrigEle, ele2->eta(), ele2->pt());
      gt.sf_eleTrig = 1 - (1-eff1)*(1-eff2);
    }
  } // done with ele trig SF
  if (gt.nLooseMuon>1 && analysis.complicatedLeptons) {
    gt.sf_muTrig = utils.getCorr(cTrigDoubleMuLeg1, gt.muonEta[0], gt.muonPt[0]) *
                   utils.getCorr(cTrigDoubleMuLeg2, gt.muonEta[1], gt.muonPt[1]);
  } else if (gt.nLooseMuon>0) {
    Muon *mu1=nullptr, *mu2=nullptr;
    if (gt.nLooseLep>0) mu1 = dynamic_cast<Muon*>(lep0);
    if (gt.nLooseLep>1) mu2 = dynamic_cast<Muon*>(lep1);
    float eff1=0, eff2=0;
    if (mu1!=nullptr && mu1->tight) {
      eff1 = utils.getCorr(
        cTrigMu,
        fabs(mu1->eta()),
        TMath::Max((float)26.,TMath::Min((float)499.99,(float)mu1->pt()))
      );
      if (mu2!=nullptr && mu2->tight)
        eff2 = utils.getCorr(
          cTrigMu,
          fabs(mu2->eta()),
          TMath::Max((float)26.,TMath::Min((float)499.99,(float)mu2->pt()))
        );
      gt.sf_muTrig = 1 - (1-eff1)*(1-eff2);
    }
  } // done with mu trig SF

  if (gt.nLoosePhoton>0 && gt.loosePho1IsTight)
    gt.sf_phoTrig = utils.getCorr(cTrigPho,gt.loosePho1Pt);

  if (analysis.vbf) {
    gt.sf_metTrigVBF = utils.getCorr(cVBF_TrigMET,gt.barrelHTMiss);
    gt.sf_metTrigZmmVBF = utils.getCorr(cVBF_TrigMETZmm,gt.barrelHTMiss);
  }
}
