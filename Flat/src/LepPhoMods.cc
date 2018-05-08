#include "../interface/LepPhoMods.h"
#include "PandaAnalysis/Utilities/interface/Helicity.h"

using namespace pa;
using namespace std;
using namespace panda;

void SimpleLeptonMod::scaleFactors()
{
  if (analysis.isData)
    return; 
  // For the hadronic analyses, store a single branch for lepton ID,
  // isolation, and tracking, computed based on the 2 leading loose
  // leptons in the event. Cool guyz get per-leg scalefactors
  // computed in ComplicatedLeptonMod
  for (int iL=0; iL!=TMath::Min(gt.nLooseLep,2); ++iL) {
    auto* lep = looseLeps.at(iL);
    float pt = lep->pt(), eta = lep->eta(), aeta = TMath::Abs(eta);
    Muon* mu = dynamic_cast<Muon*>(lep);
    if (mu!=nullptr) {
      bool isTight = mu->tight;
      if (isTight) {
        gt.sf_lepID *= utils.getCorr(cMuTightID,aeta,pt);
        gt.sf_lepIso *= utils.getCorr(cMuTightIso,aeta,pt);
      } else {
        gt.sf_lepID *= utils.getCorr(cMuLooseID,aeta,pt);
        gt.sf_lepIso *= utils.getCorr(cMuLooseIso,aeta,pt);
      }
      gt.sf_lepTrack *= utils.getCorr(cMuReco,gt.npv);
    } else {
      Electron* ele = dynamic_cast<Electron*>(lep);
      bool isTight = ele->tight;
      if (isTight) {
        gt.sf_lepID *= utils.getCorr(cEleTight,eta,pt);
      } else {
        gt.sf_lepID *= utils.getCorr(cEleVeto,eta,pt);
      }
      gt.sf_lepTrack *= utils.getCorr(cEleReco,eta,pt);
    }
  }
}

void SimpleLeptonMod::do_execute() 
{
  //electrons
  for (auto& ele : event.electrons) {
    float pt = ele.pt(); float eta = ele.eta(); float aeta = fabs(eta);
    if (!ele.veto) 
      continue; 
    if (pt<10 || aeta>2.5) 
      continue;
    if (!ElectronIP(ele.eta(),ele.dxy,ele.dz)) 
      continue;
    int iL = gt.nLooseElectron;
    bool isFake   = ele.hltsafe;
    bool isMedium = ele.medium;
    bool isTight  = ele.tight && pt>40 && aeta<2.5;
    bool isDxyz   = true; // already selected on this 
    if (isTight) gt.nTightElectron++;
    int eleSelBit            = kLoose;
    if (isFake  ) eleSelBit |= kFake;
    if (isMedium) eleSelBit |= kMedium;
    if (isTight ) eleSelBit |= kTight;
    if (isDxyz  ) eleSelBit |= kDxyz;
    if (ele.mvaWP90) eleSelBit |= kEleMvaWP90;
    if (ele.mvaWP80) eleSelBit |= kEleMvaWP80;
    gt.electronPt[iL]           = pt;
    gt.electronEta[iL]          = eta;
    gt.electronPhi[iL]          = ele.phi();
    gt.electronSelBit[iL]       = eleSelBit;
    gt.electronPdgId[iL]        = ele.charge*-11;
    looseLeps.push_back(&ele);
    matchLeps.push_back(&ele); 
    gt.nLooseElectron++;
    if (gt.nLooseElectron>=2) 
      break;
  }
  // muons
  for (auto& mu : event.muons) {
    float pt = mu.pt(); float eta = mu.eta(); float aeta = fabs(eta);
    if (!mu.loose) 
      continue; // loose ID 
    if (pt<10 || aeta>2.4) 
      continue;
    if (mu.combIso() > 0.25*pt)
      continue; // loose iso
    bool isFake   = mu.tight  && (mu.combIso() < 0.4*pt) && (mu.chIso < 0.4*pt); // what the hell is this ID?
    bool isMedium = mu.medium && (mu.combIso() < 0.15*pt);
    bool isTight  = mu.tight  && (mu.combIso() < 0.15*pt) && pt>20 && aeta<2.4;
    bool isDxyz   = MuonIP(mu.dxy,mu.dz);
    if (isTight) gt.nTightMuon++;
    int muSelBit            = kLoose;
    if (isFake  ) muSelBit |= kFake;
    if (isMedium) muSelBit |= kMedium;
    if (isTight ) muSelBit |= kTight;
    if (isDxyz  ) muSelBit |= kDxyz;
    int iL=gt.nLooseMuon;
    gt.muonPt[iL]                   = pt;
    gt.muonEta[iL]                  = eta;
    gt.muonPhi[iL]                  = mu.phi();
    gt.muonSelBit[iL]               = muSelBit;
    gt.muonPdgId[iL]                = mu.charge*-13;
    looseLeps.push_back(&mu);
    matchLeps.push_back(&mu); 
    TVector2 vMu; vMu.SetMagPhi(pt,mu.phi());
    METLOOP {
      (*jesShifts)[shift].vpfMETNoMu += vMu;
    }
    gt.nLooseMuon++;
    if (gt.nLooseMuon>=2) 
      break;
  }
  METLOOP {
    gt.pfmetnomu[shift] = (*jesShifts)[shift].vpfMETNoMu.Mod();
  }

  // now consider all leptons
  gt.nLooseLep = looseLeps.size();
  gt.nTightLep = gt.nTightElectron + gt.nTightMuon;

  if (gt.nLooseLep>0) {
    auto ptsort([](Lepton const* l1, Lepton const* l2)->bool {
      return l1->pt() > l2->pt();
    });
    int nToSort = TMath::Min(4,gt.nLooseLep);
    std::partial_sort(looseLeps.begin(),looseLeps.begin()+nToSort,looseLeps.end(),ptsort);

    Lepton* lep1 = looseLeps[0];
    METLOOP {
      gt.mT[shift] = MT(lep1->pt(),lep1->phi(),gt.pfmet[shift],gt.pfmetphi[shift]);
    }
  }  

  for (int i = 0; i != min(4, gt.nLooseLep); ++i) {
    Muon *mu = dynamic_cast<Muon*>(looseLeps[i]);
    if (mu != nullptr) {
      lepPdgId[i] = mu->charge * -13;
    } else {
      Electron *ele = dynamic_cast<Electron*>(looseLeps[i]);
      lepPdgId[i] = ele->charge * -11;
    }
  }
  
  if (gt.nLooseLep>1 && lepPdgId[0]+lepPdgId[1]==0) {
    TLorentzVector v1,v2;
    Lepton *lep1=looseLeps[0], *lep2=looseLeps[1];
    v1.SetPtEtaPhiM(lep1->pt(),lep1->eta(),lep1->phi(),lep1->m());
    v2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());
    dilep = v1+v2;
    gt.diLepMass = dilep.M();
  } else {
    gt.diLepMass = -1;
  }

  scaleFactors();
}

void ComplicatedLeptonMod::do_readData(TString dirPath)
{
  rochesterCorrection = new RoccoR(Form("%s/rcdata.2016.v3", dirPath.Data()));
}

void ComplicatedLeptonMod::do_execute()
{
  for (auto& ele : event.electrons) {
    float pt = ele.pt(); float eta = ele.eta(); float aeta = fabs(eta);
    if (analysis.hbb) {
      // Use the unsmeared/uncorrected electron pT for this loose isolation cut
      // because the electron pT assignment can be funky if it's inside a jet
      if (pt<7 || aeta>2.4 || fabs(ele.dxy)>0.05 || fabs(ele.dz)>0.2 || ele.combIso()>0.4*pt) 
        continue;
    } else {
      if (pt<10 || aeta>2.5 || !ele.veto) 
        continue;
    }
    pt = ele.smearedPt;
    ele.setPtEtaPhiM(pt,eta,ele.phi(),511e-6);
    matchLeps.push_back(&ele);
    bool isFake   = ele.hltsafe;
    bool isMedium = ele.medium;
    bool isTight  = ele.tight;
    bool isDxyz   = ElectronIP(ele.eta(),ele.dxy,ele.dz);
    bool eleMVAPresel = 
      pt > 15 && ((
        aeta < 1.4442 && 
        ele.sieie < 0.012 && 
        ele.hOverE < 0.09 &&
        ele.ecalIso < 0.4*pt && 
        ele.hcalIso < 0.25*pt && 
        ele.trackIso < 0.18*pt &&
        fabs(ele.dEtaInSeed) < 0.0095 && 
        fabs(ele.dPhiIn) < 0.065
      ) || (
        aeta > 1.5660 && 
        ele.sieie < 0.033 && 
        ele.hOverE < 0.09 &&
        ele.ecalIso < 0.45*pt && 
        ele.hcalIso < 0.28*pt &&
        ele.trackIso < 0.18*pt
    ));
    int eleSelBit            = kLoose;
    if (isFake  ) eleSelBit |= kFake;
    if (isMedium) eleSelBit |= kMedium;
    if (isTight ) eleSelBit |= kTight;
    if (isDxyz  ) eleSelBit |= kDxyz;
    if (ele.mvaWP90 && eleMVAPresel) eleSelBit |= kEleMvaWP90;
    if (ele.mvaWP80 && eleMVAPresel) eleSelBit |= kEleMvaWP80;

    // For HBB analysis, need to count as loose leptons, only the MVA90 ones
    if (analysis.hbb && !(ele.mvaWP90 && eleMVAPresel && pt>15))
      continue;
    
    if (isTight) gt.nTightElectron++;
    int iL=gt.nLooseElectron;
    gt.electronPt[iL]           = pt;
    gt.electronEta[iL]          = eta;
    gt.electronPhi[iL]          = ele.phi();
    gt.electronD0[iL]           = ele.dxy;
    gt.electronDZ[iL]           = ele.dz;
    gt.electronSfLoose[iL]      = utils.getCorr(cEleLoose, eta, pt);
    gt.electronSfMedium[iL]     = utils.getCorr(cEleMedium, eta, pt);
    gt.electronSfTight[iL]      = utils.getCorr(cEleTight, eta, pt);
    gt.electronSfMvaWP90[iL]      = utils.getCorr(cEleMvaWP90, eta, pt);
    gt.electronSfMvaWP80[iL]      = utils.getCorr(cEleMvaWP80, eta, pt);
    gt.electronSfUnc[iL]        = utils.getError(cEleMedium, eta, pt);
    gt.electronSfReco[iL]       = utils.getCorr(cEleReco, eta, pt);
    gt.electronSelBit[iL]       = eleSelBit;
    gt.electronPdgId[iL]        = ele.charge*-11;
    gt.electronNMissingHits[iL] = ele.nMissingHits;
    gt.electronTripleCharge[iL] = ele.tripleCharge;
    gt.electronCombIso[iL] = ele.combIso();
    looseLeps.push_back(&ele);
    // WARNING: The definition of "loose" here may not match your analysis 
    // definition of a loose electron for lepton multiplicity or jet cleaning 
    // considerations. It is the user's responsibility to make sure they are 
    // cutting on the correct multiplicity. Enough information is provided to 
    // do this downstream.
    gt.nLooseElectron++; 
    if (gt.nLooseElectron>=NLEP) 
      break;
  }

  // muons
  int rocRNGIdx = 0;
  for (auto& mu : event.muons) {
    float pt = mu.pt(); float eta = mu.eta(); float aeta = fabs(eta);
    if (pt<2 || aeta>2.4) continue;
    double ptCorrection=1;
    if (analysis.isData) { // perform the rochester correction on the actual particle
      ptCorrection=rochesterCorrection->kScaleDT((int)mu.charge, pt, eta, mu.phi(), 0, 0);
    } else if (pt>0) { // perform the rochester correction to the simulated particle
      // attempt gen-matching to a final state muon
      bool muonIsTruthMatched=false; TLorentzVector genP4; GenParticle genParticle;
      for (int iG = 0; iG != (int)genP->size() && !muonIsTruthMatched; ++iG) {
        genParticle = pToGRef((*genP)[iG]);
        if (genParticle.finalState != 1) continue;
        if (genParticle.pdgid != ((int)mu.charge) * -13) continue;
        genP4.SetPtEtaPhiM(genParticle.pt(), genParticle.eta(), genParticle.phi(), 0.106);
        double dR = genP4.DeltaR(mu.p4());
        if (dR < 0.3) muonIsTruthMatched=true;
      } if (muonIsTruthMatched) { // correct using the gen-particle pt
        double random1 = event.rng.uniform(rocRNGIdx);
        ptCorrection=rochesterCorrection->kScaleFromGenMC((int)mu.charge, 
                                                          pt, eta, mu.phi(), 
                                                          mu.trkLayersWithMmt, 
                                                          genParticle.pt(), 
                                                          random1, 0, 0);
      } else { // if gen match not found, correct the other way
        double random1 = event.rng.uniform(rocRNGIdx); 
        double random2 = event.rng.uniform(rocRNGIdx);
        ptCorrection=rochesterCorrection->kScaleAndSmearMC((int)mu.charge, 
                                                           pt, eta, mu.phi(), 
                                                           mu.trkLayersWithMmt, 
                                                           random1, random2, 0, 0);
      }
    }
    pt *= ptCorrection;
    if (analysis.hbb) {
      if (pt<5 || aeta>2.4 || !mu.loose || fabs(mu.dxy)>0.5 || fabs(mu.dz)>1.0 || mu.combIso()>0.4*pt) continue;
    } else {
      if (pt<10 || aeta>2.4 || !mu.loose) continue;
    }
    mu.setPtEtaPhiM(pt,eta,mu.phi(),0.106);
    matchLeps.push_back(&mu); // Muon jet cleaning and loose cuts are equivalent in HBB land

    bool isFake   = mu.tight  && mu.combIso()/mu.pt() < 0.4 && mu.chIso/mu.pt() < 0.4;
    bool isMedium = mu.medium && mu.combIso()/mu.pt() < 0.15;
    bool isTight  = mu.tight  && mu.combIso()/mu.pt() < 0.15;
    bool isDxyz   = MuonIP(mu.dxy,mu.dz);
    if (isTight) gt.nTightMuon++;
    int muSelBit            = kLoose;
    if (isFake  ) muSelBit |= kFake;
    if (isMedium) muSelBit |= kMedium;
    if (isTight ) muSelBit |= kTight;
    if (isDxyz  ) muSelBit |= kDxyz;
    int iL=gt.nLooseMuon;
    gt.muonPt[iL]                   = pt;
    gt.muonEta[iL]                  = eta;
    gt.muonPhi[iL]                  = mu.phi();
    gt.muonD0[iL]                   = mu.dxy;
    gt.muonDZ[iL]                   = mu.dz;
    gt.muonSfLoose[iL]              = utils.getCorr(cMuLooseID, TMath::Abs(mu.eta()), mu.pt()) 
                                        * utils.getCorr(cMuLooseIso, TMath::Abs(mu.eta()), 
                                                    mu.pt());
    gt.muonSfMedium[iL]             = utils.getCorr(cMuMediumID, TMath::Abs(mu.eta()), mu.pt());
    gt.muonSfTight[iL]              = utils.getCorr(cMuTightID, TMath::Abs(mu.eta()), mu.pt())  
                                        * utils.getCorr(cMuTightIso, TMath::Abs(mu.eta()), mu.pt());
    gt.muonSfUnc[iL]                = utils.getError(cMuMediumID , TMath::Abs(mu.eta()), mu.pt());
    gt.muonSfReco[iL]               = utils.getCorr(cMuReco, mu.eta());
    gt.muonSelBit[iL]               = muSelBit;
    gt.muonPdgId[iL]                = mu.charge*-13;
    gt.muonIsSoftMuon[iL]           = mu.soft;
    gt.muonCombIso[iL]              = mu.combIso();
    looseLeps.push_back(&mu);
    TVector2 vMu; vMu.SetMagPhi(pt,mu.phi());
    METLOOP {
      (*jesShifts)[shift].vpfMETNoMu += vMu;
    }
    // WARNING: The definition of "loose" here may not match your analysis 
    // definition of a loose muon for lepton multiplicity or jet cleaning 
    // considerations. It is the user's responsibility to make sure they are 
    // cutting on the correct multiplicity. Enough information is provided to 
    // do this downstream.
    gt.nLooseMuon++;
    if (gt.nLooseMuon>=NLEP) 
      break;
  }
  METLOOP {
    gt.pfmetnomu[shift] = (*jesShifts)[shift].vpfMETNoMu.Mod();
  }

  // now consider all leptons
  gt.nLooseLep = looseLeps.size();
  gt.nTightLep = gt.nTightElectron + gt.nTightMuon;
  if (gt.nLooseLep>0) {
    auto ptsort([](Lepton const* l1, Lepton const* l2)->bool {
      return l1->pt() > l2->pt();
    });
    int nToSort = TMath::Min(12,gt.nLooseLep);
    std::partial_sort(looseLeps.begin(),looseLeps.begin()+nToSort,looseLeps.end(),ptsort);

    Lepton* lep1 = looseLeps[0];
    METLOOP {
      gt.mT[shift] = MT(lep1->pt(),lep1->phi(),gt.pfmet[shift],gt.pfmetphi[shift]);
    }
  }

  for (int i = 0; i != min(4, gt.nLooseLep); ++i) {
    Muon *mu = dynamic_cast<Muon*>(looseLeps[i]);
    if (mu != nullptr) {
      lepPdgId[i] = mu->charge * -13;
    } else {
      Electron *ele = dynamic_cast<Electron*>(looseLeps[i]);
      lepPdgId[i] = ele->charge * -11;
    }
  }
  
  if (gt.nLooseLep>1) {
    TLorentzVector v1,v2;
    Lepton *lep1=looseLeps[0], *lep2=looseLeps[1];
    v1.SetPtEtaPhiM(lep1->pt(),lep1->eta(),lep1->phi(),lep1->m());
    v2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());
    dilep = v1+v2;
    if (lepPdgId[0]+lepPdgId[1]==0)  // Strict sign/flavor 
      gt.diLepMass = dilep.M();
  
    // Also allow for opposite flavor or same sign selection:
    gt.ZBosonPt  = dilep.Pt();
    gt.ZBosonEta = dilep.Eta();
    gt.ZBosonPhi = dilep.Phi();
    gt.ZBosonM   = dilep.M();
    gt.ZBosonLep1CosThetaCS = CosThetaCollinsSoper(v1,v2);
    // ZBosonLep1CosThetaStar, ZBosonLep1CosThetaStarFJ
    // are not calculated here, we need the Z(ll)H(bb) system in JetsMods
    
  } 

}


void InclusiveLeptonMod::do_execute()
{
  // "Inclusive very loose selection"
  // https://github.com/vhbb/cmssw/blob/vhbbHeppy80X/PhysicsTools/Heppy/python/analyzers/objects/LeptonAnalyzer.py
  for (auto& ele : event.electrons) {
    float pt = ele.pt(); float eta = ele.eta(); float aeta = fabs(eta);
    if (pt<5 || aeta>2.5) 
      continue;
    if (fabs(ele.dxy)>0.5 || fabs(ele.dz)>1.0) 
      continue;
    if (ele.nMissingHits>1) 
      continue;
    inclusiveLeps.push_back(&ele);
  }
  for (auto& mu : event.muons) {
    if (!mu.loose) 
      continue;
    float pt = mu.pt(); float eta = mu.eta(); float aeta = fabs(eta);
    if (pt<3 || aeta>2.4) 
      continue;
    if (fabs(mu.dxy)>0.5 || fabs(mu.dz)>1.0) 
      continue;
    inclusiveLeps.push_back(&mu);
  }
}

void SimplePhotonMod::scaleFactors()
{
  if (analysis.isData)
    return; 
  if (gt.nLoosePhoton < 1)
    return;
  float pt = gt.loosePho1Pt, eta = gt.loosePho1Eta;
  if (gt.loosePho1IsTight)
    gt.sf_pho = utils.getCorr(cPho,eta,pt);
  if (analysis.isData && pt>175) {
    gt.sf_phoPurity = utils.getCorr(cPhoFake, pt);
  }
}

void SimplePhotonMod::do_execute()
{
  for (auto& pho : event.photons) {
    if (!pho.loose || !pho.csafeVeto)
      continue;
    float pt = pho.pt();
    if (pt<1) 
      continue;
    float eta = pho.eta(), phi = pho.phi();
    if (pt<15 || fabs(eta)>2.5)
      continue;
    loosePhos.push_back(&pho);
    gt.nLoosePhoton++;
    if (gt.nLoosePhoton==1) {
      gt.loosePho1Pt = pt;
      gt.loosePho1Eta = eta;
      gt.loosePho1Phi = phi;
    }
    if ( pho.medium &&
         pt>175 ) { // apply eta cut offline
      if (gt.nLoosePhoton==1)
        gt.loosePho1IsTight=1;
      gt.nTightPhoton++;
      tightPhos.push_back(&pho);
    }
  }

  scaleFactors();
}


void ComplicatedPhotonMod::do_execute()
{
  for (auto& pho : event.photons) {
    if (!pho.medium)
      continue;
    float pt = pho.pt();
    if (pt<1) 
      continue;
    float eta = pho.eta(), phi = pho.phi();
    if (pt<25 || fabs(eta)>2.5)
      continue;
    if (isMatched(matchLeps,0.16,pho.eta(),pho.phi()))
      continue;
    loosePhos.push_back(&pho);
    gt.nLoosePhoton++;
    if (gt.loosePho1Pt < pt) {
      gt.loosePho1Pt = pt;
      gt.loosePho1Eta = eta;
      gt.loosePho1Phi = phi;
      int phoSelBit = 0;
      // this is always true as of now, but safer to have it like this
      if (pho.medium)                 phoSelBit |= pMedium; 
      if (pho.tight)                  phoSelBit |= pTight;
      if (pho.highpt)                 phoSelBit |= pHighPt;
      if (pho.csafeVeto)              phoSelBit |= pCsafeVeto;
      if (pho.pixelVeto)              phoSelBit |= pPixelVeto;
      if (!pfChargedPhotonMatch(pho)) phoSelBit |= pTrkVeto;
      gt.loosePho1SelBit = phoSelBit;
      if (pho.medium && pho.csafeVeto && pho.pixelVeto) gt.loosePho1IsTight = 1;
      else                                              gt.loosePho1IsTight = 0;
    }
    if ( pho.medium && pho.csafeVeto && pho.pixelVeto) { // apply eta cut offline
      gt.nTightPhoton++;
      tightPhos.push_back(&pho);
    }
  }

  scaleFactors();
}

bool ComplicatedPhotonMod::pfChargedPhotonMatch(const Photon& photon)
{
  double matchedRelPt = -1.;

  for (auto& cand : event.pfCandidates) {
    if (cand.q() == 0) 
      continue;

    double dr(cand.dR(photon));
    double rawPt = photon.pt();
    double relPt(cand.pt() / rawPt);
    if (dr < 0.1 && relPt > matchedRelPt) {
      matchedRelPt = relPt;
    }

  }

  return (matchedRelPt > 0.6);
}

void TauMod::do_execute()
{
  for (auto& tau : event.taus) {
    if (analysis.vbf) {
      if (!tau.decayMode || !tau.decayModeNew)
        continue;
      if (!tau.looseIsoMVAOld)
        continue;
    } else {
      if (!tau.decayMode || !tau.decayModeNew)
        continue;
      if (!tau.looseIsoMVA)
        continue;
    }
    if (tau.pt()<18 || fabs(tau.eta())>2.3)
      continue;
    if (isMatched(matchLeps,0.16,tau.eta(),tau.phi()))
      continue;
    gt.nTau++;
  }
}

void GenLepMod::do_execute()
{
  gt.genTauPt = -1;
  gt.genElectronPt = -1;
  gt.genMuonPt = -1;
  GenParticle *tau = NULL;
  bool foundTauLeptonic = false; 
  for (auto* genptr : *genP) {
    auto& gen = pToGRef(genptr);
    int apdgid = abs(gen.pdgid);
    float pt = gen.pt();
    bool isEmu = false; 

    if (apdgid == 11 && pt > gt.genElectronPt) {
      gt.genElectronPt = pt; 
      gt.genElectronEta = gen.eta(); 
      isEmu = true; 
    }
    
    if (apdgid == 13 && pt > gt.genMuonPt) {
      gt.genMuonPt = pt; 
      gt.genMuonEta = gen.eta(); 
      isEmu = true; 
    }

    if (isEmu && !foundTauLeptonic && tau) {
      const GenParticle *parent = &gen;
      while (parent->parent.isValid()) {
        parent = parent->parent.get();
        if (parent == tau) {
          foundTauLeptonic = true; 
          gt.genTauPt = -1; 
          gt.genTauEta = -1;
          break;
        }
      }
    }

    if (!foundTauLeptonic && apdgid == 15 && pt > gt.genTauPt
        && ((gen.statusFlags & (1 << GenParticle::kIsHardProcess)) != 0 
            || (gen.statusFlags & (1 << GenParticle::kFromHardProcessBeforeFSR)) != 0 
            || ((gen.statusFlags & (1 << GenParticle::kIsDecayedLeptonHadron)) != 0 
                && (gen.statusFlags & (1 << GenParticle::kFromHardProcess)) != 0
                )
            )
        ) 
    {
      gt.genTauPt = pt; 
      gt.genTauEta = gen.eta();
    }
  }
}
