#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>
#include "PandaAnalysis/Utilities/interface/NeutrinoSolver.h"
#include "PandaAnalysis/Utilities/interface/Helicity.h"
#define TWOPI 6.28318531

using namespace panda;
using namespace std;

float centralOnly(float x, float aeta, float def = -1) 
{
  return  aeta < 2.4 ? x : -1;
}

bool csvLoose(float csv) 
{
  return csv > 0.5426;
}

bool csvMed(float csv) 
{
  return csv > 0.8484;
}

PandaAnalyzer::JetWrapper shiftJet(const panda::Jet& jet, shiftjes shift) 
{
  float pt;
  switch (shift) {
    case shiftjes::kNominal:
      pt = jet.pt();
      break;
    case shiftjes::kJESUp:
      pt = jet.ptCorrUp;
      break;
    case shiftjes::kJESDown:
      pt = jet.ptCorrDown;
      break;
    default:
      PError("shiftJet", "Unknown JES type!");
      exit(1);
  }
  return PandaAnalyzer::JetWrapper(pt, jet);
}

// set up JES readers
// Responsible: S. Narayanan
void PandaAnalyzer::SetupJES()
{
  if (uncReader==0) {
    if (isData) {
      TString thisEra = eras.getEra(gt->runNumber);
      for (auto &iter : ak8UncReader) {
        if (! iter.first.Contains("data"))
          continue;
        if (iter.first.Contains(thisEra)) {
          uncReader = iter.second;
          uncReaderAK4 = ak4UncReader[iter.first];
          scaleReaderAK4 = ak4ScaleReader[iter.first];
          break;
        }
      }
    } else {
      uncReader = ak8UncReader["MC"];
      uncReaderAK4 = ak4UncReader["MC"];
      scaleReaderAK4 = ak4ScaleReader["MC"];
    }
  }
}

void PandaAnalyzer::jetPartonFlavor(const Jet& jet, int& flavor, float& genpt) {
  // here we try to match to hard partons
  flavor=0; genpt=0; // defaults
  for (auto* genptr : validGenP) {
    auto& gen = pToGRef(genptr);
    int apdgid = abs(gen.pdgid);
    if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
      continue;
    double dr2 = DeltaR2(jet.eta(),jet.phi(),gen.eta(),gen.phi());
    if (dr2<0.09) {
      genpt = gen.pt();
      if (apdgid==4 || apdgid==5) {
        flavor=apdgid;
        return;
      } else {
        flavor=0;
      }
    }
  } 
}

void PandaAnalyzer::jetClusteredFlavor(const Jet& jet, int& flavor, float& genpt) {
  // or we can match to gen jets (probably better)
  flavor=0; genpt=0; // defaults
  for (auto &gen : event.ak4GenJets) {
    if (DeltaR2(gen.eta(), gen.phi(), jet.eta(), jet.phi()) < 0.09) {
      int apdgid = abs(gen.pdgid);
      genpt = gen.pt();
      if (apdgid == 4 || apdgid == 5) {
        flavor = apdgid;
        return;
      } else {
        flavor = 0;
      }
    }
  }
}

// basic jet information
// Responsible: S. Narayanan
void PandaAnalyzer::JetBasics() 
{
  float maxJetEta = analysis->vbf ? 4.7 : 4.5;
  int nJetDPhi = analysis->vbf ? 4 : 5;
  float minMinJetPt = min(minJetPt, minBJetPt);

  int maxshift = analysis->varyJES ? jes2i(shiftjes::N) : 1;
  TLorentzVector vBarrelJets;

  for (int shift = jes2i(shiftjes::kNominal); shift != maxshift; ++shift) {
    JESHandler& jets = jesShifts[shift];
    for (auto& jw : jets.all) {
      auto& jet = jw.get_base();
      float aeta = abs(jet.eta());
      float pt = jw.pt;
      if (aeta > maxJetEta || pt < minMinJetPt) 
        continue;
      if (IsMatched(&matchLeps,0.16,jet.eta(),jet.phi()))
        continue;
      if (!analysis->hbb && IsMatched(&matchPhos,0.16,jet.eta(),jet.phi()))
        continue;
      if ((analysis->vbf || analysis->hbb) && !jet.loose)
        continue;

      if (analysis->jetFlavorPartons)
        jetPartonFlavor(jet, jw.flavor, jw.genpt);
      else if (analysis->jetFlavorJets)
        jetClusteredFlavor(jet, jw.flavor, jw.genpt);

      float csv = centralOnly(jet.csv, aeta);
      float cmva = centralOnly(jet.cmva, aeta);

      if (pt > minBJetPt && aeta < 2.4) { // b jets
        if (csvLoose(csv)) { 
          ++(gt->jetNBtags[shift]);
          if (csvMed(csv)) {
            ++(gt->jetNMBtags[shift]);
          }
        }

        jets.bcand.push_back(&jw); 
      } 


      if (pt > minJetPt) {
        if ((analysis->hbb || analysis->monoh) && jets.cleaned.size() >= NJET)
          continue;

        jets.cleaned.push_back(&jw);
        
        if (jets.cleaned.size() < 3) {      
          if (GetCorr(cBadECALJets, jet.eta(), jet.phi()) > 0)
            gt->badECALFilter = 0;
        }

        if (analysis->fatjet) {
          IsoJet(jw, jets);
          if (jw.iso && csvLoose(csv))
            ++(gt->isojetNBtags[shift]);
        }

        if (aeta < 2.4) {
          jets.central.push_back(&jw);
          
          int njet = jets.central.size();
          gt->nJet[shift] = njet;
          if (shift == jes2i(shiftjes::kNominal)) {
            if (njet < 2) {
              gt->jetPt[njet] = pt;
              gt->jetEta[njet] = jet.eta();
              gt->jetPhi[njet] = jet.phi();
              gt->jetCSV[njet] = csv;
              gt->jetIsTight[njet] = jet.monojet ? 1 : 0;
              gt->jetFlav[njet] = jw.flavor;
              gt->jetIsIso[njet] = jw.iso ? 1 : 0;
            }
          }
        }

        if (shift == jes2i(shiftjes::kNominal) && aeta < 3.0) {
          gt->barrelHT += jet.pt();
          vBarrelJets += jet.p4();
        }

        int njet = jets.cleaned.size();
        if (njet < 2 || ((analysis->hbb || analysis->monoh) && njet < NJET)) {
          gt->jotPt[njet][shift] = pt;
          gt->jotEta[njet][shift] = jet.eta();
          gt->jotPhi[njet][shift] = jet.phi();
          gt->jotE[njet][shift] = jet.e() * jw.scale();
          gt->jotCSV[njet][shift] = csv;
          gt->jotCMVA[njet][shift] = cmva;
          gt->jotVBFID[njet][shift] = (aeta < 2.4) ? (jet.monojet ? 1 : 0) : 1;

          if (analysis->bjetRegression) {
            JetBRegressionInfo(jet, njet, shift);
          }
        }

        TLorentzVector vJet; vJet.SetPtEtaPhiM(pt, jet.eta(), jet.phi(), jet.m());
        if (njet < nJetDPhi) {
          gt->dphipuppimet[shift] = min(fabs(vJet.DeltaPhi(jets.vpuppiMET)), (double)gt->dphipuppimet[shift]);
          gt->dphipfmet[shift]    = min(fabs(vJet.DeltaPhi(jets.vpfMET)),    (double)gt->dphipfmet[shift]);
          if (analysis->recoil) {
            gt->dphipuppiUA[shift] = min(fabs(vJet.DeltaPhi(jets.vpuppiUA)), (double)gt->dphipuppiUA[shift]);
            gt->dphipuppiUW[shift] = min(fabs(vJet.DeltaPhi(jets.vpuppiUW)), (double)gt->dphipuppiUW[shift]);
            gt->dphipuppiUZ[shift] = min(fabs(vJet.DeltaPhi(jets.vpuppiUZ)), (double)gt->dphipuppiUZ[shift]);
            gt->dphipfUA[shift] = min(fabs(vJet.DeltaPhi(jets.vpfUA)), (double)gt->dphipfUA[shift]);
            gt->dphipfUW[shift] = min(fabs(vJet.DeltaPhi(jets.vpfUW)), (double)gt->dphipfUW[shift]);
            gt->dphipfUZ[shift] = min(fabs(vJet.DeltaPhi(jets.vpfUZ)), (double)gt->dphipfUZ[shift]);
          }
        }
      }
    }

    gt->nJot[shift] = jets.cleaned.size();
    switch (gt->whichRecoil) {
      case 0: // MET
        gt->dphipuppiU[shift] = gt->dphipuppimet[shift];
        gt->dphipfU[shift] = gt->dphipfmet[shift];
        break;
      case -1: // photon
        gt->dphipuppiU[shift] = gt->dphipuppiUA[shift];
        gt->dphipfU[shift] = gt->dphipfUA[shift];
        break;
      case 1:
        gt->dphipuppiU[shift] = gt->dphipuppiUW[shift];
        gt->dphipfU[shift] = gt->dphipfUW[shift];
        break;
      case 2:
        gt->dphipuppiU[shift] = gt->dphipuppiUZ[shift];
        gt->dphipfU[shift] = gt->dphipfUZ[shift];
        break;
      default: // c'est impossible !
        break;
    }

    if (analysis->vbf) 
      JetVBFSystem(shift);

    if (analysis->monoh || analysis->hbb)
      JetHbbReco(shift);

  } // shift loop 

  gt->barrelHTMiss = vBarrelJets.Pt();

  tr->TriggerEvent("jets");

}

// stuff for b-jet energy regression
// Responsible: B. Maier, D. Hsu
void PandaAnalyzer::JetBRegressionInfo(const panda::Jet& jet, int N, int shift)
{
  gt->jotEMF[N][shift] = jet.cef + jet.nef;
  gt->jotHF[N][shift] = jet.chf + jet.nhf;
  gt->jotLep1Pt[N][shift] = 0;
  gt->jotTrk1Pt[N][shift] = 0;
  gt->jotNLep[N][shift] = 0;
  for (const panda::ConstRef<panda::PFCand> &c_iter : jet.constituents) {
    if (!c_iter.isValid())
      continue;
    auto *pf = c_iter.get();
    if (pf->q() != 0) {
      float pt = pf->pt();
      gt->jotTrk1Pt[N][shift] = max(pt, gt->jotTrk1Pt[N][shift]);
      int pdgid = abs(pf->pdgId());
      if (pdgid == 11 || pdgid == 13) {
        gt->jotNLep[N][shift]++;
        if (pt > gt->jotLep1Pt[N][shift]) {
          gt->jotLep1Pt[N][shift] = pt;
          gt->jotLep1PtRel[N][shift] = pf->p4().Perp(jet.p4().Vect());
          gt->jotLep1DeltaR[N][shift] = sqrt(DeltaR2(pf->eta(), pf->phi(), jet.eta(), jet.phi()));
        }
      }
    }
  }

  auto& vert = jet.secondaryVertex;
  if (vert.isValid()) {
    gt->jotVtxPt[N][shift] = vert->pt();
    gt->jotVtxMass[N][shift] = vert->m();
    gt->jotVtx3DVal[N][shift] = vert->vtx3DVal;
    gt->jotVtx3DErr[N][shift] = vert->vtx3DeVal;
    gt->jotVtxNtrk[N][shift] = vert->ntrk;
  }

  tr->TriggerSubEvent("b-jet reg info");
}

// isolated jet information
// Responsible: S. Narayanan
void PandaAnalyzer::IsoJet(JetWrapper& jw, JESHandler& jets) 
{
  auto& jet = jw.get_base();
  float maxIsoEta = (analysis->monoh) ? 4.5 : 2.5;
  bool isIsoJet = ( 
        gt->nFatjet == 0 || 
        fabs(jet.eta() < maxIsoEta && 
        DeltaR2(gt->fjEta,gt->fjPhi,jet.eta(),jet.phi()) > FATJETMATCHDR2) 
      ); 

  jw.iso = isIsoJet; 

  if (isIsoJet) 
    jets.iso.push_back(&jw);

  tr->TriggerSubEvent("iso jets");
}

// vary jet energy scales for jet
// Responsible: S. Narayanan
void PandaAnalyzer::JetVaryJES()
{
  JESLOOP {
    auto& jets = jesShifts[shift];
    jets.reserve(ak4jets->size());
    for (auto &j : *ak4jets) {
      jets.all.push_back(shiftJet(j, i2jes(shift)));
    }
    if (shift != jes2i(shiftjes::kNominal)) {
      std::sort(jets.all.begin(), jets.all.end(),
                [](JetWrapper x, JetWrapper y) { return x.pt > y.pt; });
    }
  }
  tr->TriggerSubEvent("vary jet JES");
}

// variables of the entire VBF system
// Responsible: S. Narayanan
void PandaAnalyzer::JetVBFSystem(int shift) 
{
  auto& jets = jesShifts[shift];
  if (jets.cleaned.size() > 1) {
    TLorentzVector v0 = jets.cleaned[0]->p4();
    TLorentzVector v1 = jets.cleaned[1]->p4();
    gt->jot12Mass[shift] = (v0 + v1).M();
    gt->jot12DPhi[shift] = v0.DeltaPhi(v1);
    gt->jot12DEta[shift] = fabs(v0.Eta() - v1.Eta());
  }

  tr->TriggerSubEvent("VBF jet system");
}

// H->bb reconstruction
// Responsible: B. Maier, D.Hsu, S.Narayanan
void PandaAnalyzer::JetHbbReco(int shift) 
{
  auto& jets = jesShifts[shift];
  if (jets.central.size() < 2)
    return;

  vector<const JetWrapper*> btagsorted = jets.central; // copy
  sort(btagsorted.begin(), btagsorted.end(),
       analysis->useCMVA ? 
        [](const JetWrapper *x, const JetWrapper *y) { return x->base->cmva > y->base->cmva; } :   
        [](const JetWrapper *x, const JetWrapper *y) { return x->base->csv > y->base->csv; }
      );
  map<const JetWrapper*, int> order; // needed for output indexing
  for (int i = 0; i != (int)jets.cleaned.size(); ++i) 
    order[jets.cleaned[i]] = i;

  const JetWrapper* bjets[2] = {btagsorted[0], btagsorted[1]};

  gt->hbbjtidx[0][shift] = order[bjets[0]];
  gt->hbbjtidx[1][shift] = order[bjets[1]];

  TLorentzVector hbbd[2] = { bjets[0]->p4(), bjets[1]->p4() };
  TLorentzVector hbbsystem = hbbd[0] + hbbd[1];

  gt->hbbpt[shift] = hbbsystem.Pt();
  gt->hbbeta[shift] = hbbsystem.Eta();
  gt->hbbphi[shift] = hbbsystem.Phi();
  gt->hbbm[shift] = hbbsystem.M();

  tr->TriggerSubEvent("Bare Hbb reco");

  TLorentzVector hbbd_corr[2];
  if (analysis->bjetRegression && gt->hbbm[shift] > 0) {
    for (int i = 0; i<2; i++) {
      int idx = gt->hbbjtidx[i][shift];
      // Shifted values for the jet energies to perform the b-jet regression
      bjetreg_vars[0] = gt->jotPt[idx][shift];
      bjetreg_vars[1] = gt->jotEta[idx][shift];
      bjetreg_vars[2] = gt->jotTrk1Pt[idx][shift];
      bjetreg_vars[3] = gt->jotLep1Pt[idx][shift];
      bjetreg_vars[4] = gt->jotEMF[idx][shift];
      bjetreg_vars[5] = gt->jotHF[idx][shift];
      bjetreg_vars[6] = gt->jotLep1DeltaR[idx][shift];
      bjetreg_vars[7] = gt->jotLep1PtRel[idx][shift];
      bjetreg_vars[8] = gt->jotVtxPt[idx][shift];
      bjetreg_vars[9] = gt->jotVtxMass[idx][shift];
      bjetreg_vars[10]= gt->jotVtx3DVal[idx][shift];
      bjetreg_vars[11]= gt->jotVtx3DErr[idx][shift];
      bjetreg_vars[12]= gt->jotVtxNtrk[idx][shift];
      bjetreg_vars[13]= hbbd[i].Et();
      bjetreg_vars[14]= hbbd[i].Mt();

      gt->hbbjtRegFac[i][shift] = (bjetregReader->EvaluateRegression("BDT method"))[0];
      hbbd_corr[i].SetPtEtaPhiM(
            gt->hbbjtRegFac[i][shift] * gt->jotPt[idx][shift],
            gt->jotEta[idx][shift],
            gt->jotPhi[idx][shift],
            gt->jotM[idx][shift]
          );
    }

    TLorentzVector hbbsystem_corr = hbbd_corr[0] + hbbd_corr[1];
    gt->hbbm_reg[shift] = hbbsystem_corr.M(); 
    gt->hbbpt_reg[shift] = hbbsystem_corr.Pt();

    tr->TriggerSubEvent("Regr. Hbb reco");
  } // regression

  if (gt->hbbm > 0) {
    gt->hbbCosThetaJJ[shift] = hbbsystem.CosTheta();
    float csj1;
    if (analysis->bjetRegression) {
      if (hbbd_corr[0].Pt() > hbbd_corr[1].Pt())
        csj1 = CosThetaCollinsSoper(hbbd_corr[0], hbbd_corr[1]);
      else
        csj1 = CosThetaCollinsSoper(hbbd_corr[1], hbbd_corr[0]);
    } else {
      if (hbbd[0].Pt() > hbbd[1].Pt())
        csj1 = CosThetaCollinsSoper(hbbd[0], hbbd[1]);
      else
        csj1 = CosThetaCollinsSoper(hbbd[1], hbbd[0]);
    }
    gt->hbbCosThetaCSJ1[shift] = csj1;
    tr->TriggerSubEvent("Hbb spin correl.");
  }

  if (gt->hbbm > 0 && gt->nLooseLep > 0 ) {
    TLorentzVector leptonP4, metP4, nuP4, *jet0P4, *jet1P4, WP4, topP4;
    float dRJet0W, dRJet1W;
    bool jet0IsCloser;
    leptonP4=looseLeps[0]->p4();

    metP4.SetPtEtaPhiM(gt->pfmet[shift], 0, gt->pfmetphi[shift], 0);
    nuP4 = getNu4Momentum( leptonP4, metP4 );
    WP4 = leptonP4 + nuP4;

    // If using b-jet regression, use the regressed jets for the top mass reconstruction
    // Otherwise, use the un regressed jets
    if (analysis->bjetRegression) {
      jet0P4 = &hbbd_corr[0]; jet1P4 = &hbbd_corr[1];
    } else {
      jet0P4 = &hbbd[0]; jet1P4 = &hbbd[1];
    }

    dRJet0W=jet0P4->DeltaR(leptonP4); 
    dRJet1W=jet1P4->DeltaR(leptonP4); 
    jet0IsCloser = (dRJet0W < dRJet1W);

    topP4 = jet0IsCloser ? (*jet0P4)+WP4 : (*jet1P4)+WP4;
    gt->topMassLep1Met[shift] = topP4.M();
    gt->topWBosonCosThetaCS[shift] = CosThetaCollinsSoper(WP4, jet0IsCloser ? *jet0P4 : *jet1P4);
    gt->topWBosonPt  = WP4.Pt();
    gt->topWBosonEta = WP4.Eta();
    gt->topWBosonPhi = WP4.Phi();

    tr->TriggerSubEvent("Top(bW) reco");
  }

  
  tr->TriggerEvent("Hbb");
}

// clustering gen jets with neutrinos
// Responsible: S. Narayanan
void PandaAnalyzer::GenJetsNu()
{

  std::vector<fastjet::PseudoJet> finalStates;
  std::vector<const panda::GenParticle*> bcs;
  for (auto* pptr : validGenP) {
    auto& p = pToGRef(pptr);
    if (p.finalState && p.pt() > 0.001) {
      finalStates.emplace_back(p.px(), p.py(), p.pz(), p.e());
      continue;
    }
    int apdgid = abs(p.pdgid);
    if (apdgid == 4 || apdgid == 5) {
      bcs.push_back(&p);
    }
  }

  fastjet::ClusterSequenceArea seq(finalStates, *jetDefGen, *areaDef);
  std::vector<fastjet::PseudoJet> allJets(seq.inclusive_jets(0.01));

  genJetsNu.reserve(allJets.size());
  for (auto &pj : allJets) {
    genJetsNu.emplace_back(panda::GenJet());
    genJetsNu.back().setXYZE(pj.px(), pj.py(), pj.pz(), pj.e());
    int flavor = 0;
    for (auto *bc : bcs) {
      if (DeltaR2(pj.eta(), pj.phi(), bc->eta(), bc->phi()) < 0.09) {
        flavor = abs(bc->pdgid);
        break;
      }
    }
    genJetsNu.back().pdgid = flavor;
  }
  tr->TriggerEvent("gen jets nu");
}

// soft activity of event
// Responsible: D. Hsu
void PandaAnalyzer::JetHbbSoftActivity() {
  // Soft activity
  int nomidx = jes2i(shiftjes::kNominal);

  if (gt->hbbm[nomidx] > 0.) {
    gt->sumEtSoft1=0; gt->nSoft2=0; gt->nSoft5=0; gt->nSoft10=0;
    // Define the ellipse of particles to forget about
    // https://math.stackexchange.com/questions/426150/what-is-the-general-equation-of-the-ellipse-that-is-not-in-the-origin-and-rotate 
    // ((x-h)cos(A) + (y-k)sin(A))^2 /a^2 + ((x-h)sin(A) - (y-k)cos(A))^2 /b^2 <=1
    double ellipse_cosA, ellipse_sinA, ellipse_h, ellipse_k, ellipse_a, ellipse_b; {
      double ellipse_alpha;
      float phi1=gt->jetPhi[gt->hbbjtidx[0][nomidx]], phi2=gt->jetPhi[gt->hbbjtidx[1][nomidx]];
      float eta1=gt->jetEta[gt->hbbjtidx[0][nomidx]], eta2=gt->jetEta[gt->hbbjtidx[1][nomidx]];
      double phi1MinusPhi2 = phi1-phi2;
      double eta1MinusEta2 = eta1-eta2;
      double phi1MinusPhi2MPP = TVector2::Phi_mpi_pi(phi1MinusPhi2);
      ellipse_alpha = atan2( phi1MinusPhi2, eta1MinusEta2);
      // compute delta R using already computed qty's to save time
      ellipse_a = (sqrt(pow(eta1MinusEta2,2) 
                   + pow(TVector2::Phi_mpi_pi(phi1MinusPhi2),2)) + 1.)/2.; // Major axis 2*a = dR(b,b)+1
      ellipse_b = 1./2.; // Minor axis 2*b = 1
      ellipse_h = (eta1+eta2)/2.;
      ellipse_k = TVector2::Phi_mpi_pi(phi2 + phi1MinusPhi2MPP/2.);
      ellipse_cosA = cos(ellipse_alpha);
      ellipse_sinA = sin(ellipse_alpha);
      if (DEBUG > 10) {
        PDebug("PandaAnalyzer::JetHbbSoftActivity",
               Form("Calculating ellipse with (eta1,phi1)=(%.2f,%.2f), (eta2,phi2)=(%.2f,%.2f)",
                    eta1,phi1,eta2,phi2));
        PDebug("PandaAnalyzer::JetHbbSoftActivity",
               Form("Found ellipse parameters (a,b,h,k,alpha)=(%.2f,%.2f,%.2f,%.2f,%.2f)",
                    ellipse_a,ellipse_b,ellipse_h,ellipse_k,ellipse_alpha));
      }
    }

    // Find out which PF constituents to not use
    const RefVector<PFCand> &jet1Tracks = jesShifts[nomidx].cleaned[gt->hbbjtidx[0][nomidx]]->base->constituents,
                            &jet2Tracks = jesShifts[nomidx].cleaned[gt->hbbjtidx[1][nomidx]]->base->constituents;

    // Get vector of pseudo jets for clustering
    panda::PFCandCollection &allTracks = event.pfCandidates;
    std::vector<fastjet::PseudoJet> softTracksPJ;
    softTracksPJ.reserve(allTracks.size());
    panda::PFCand *softTrack=0;
    for (auto &softTrackRef : allTracks) {
      softTrack = &softTrackRef;
      // Minimum track pT threshold (300 MeV default)
      if (softTrack->pt() < minSoftTrackPt) continue;
      // High quality track flag
      if (!softTrack->track.isValid() || !softTrack->track.get()->highPurity) continue;
      // Only consider tracks with dz < 0.2 w.r.t. the primary vertex
      if (fabs(softTrack->track.get()->dz()) > 0.2) continue;
      // Track cannot be a constituent of loose leptons or the two b-jets
      bool trackIsSpokenFor=false;
      if (!trackIsSpokenFor) for (UShort_t iJetTrack=0; iJetTrack<jet1Tracks.size(); iJetTrack++) {
        if (!jet1Tracks.at(iJetTrack).isValid()) continue;
        if (softTrack==jet1Tracks.at(iJetTrack).get()) { trackIsSpokenFor=true; break; }
      }
      if (!trackIsSpokenFor) for (UShort_t iJetTrack=0; iJetTrack<jet2Tracks.size(); iJetTrack++) {
        if (!jet2Tracks.at(iJetTrack).isValid()) continue;
        if (softTrack==jet2Tracks.at(iJetTrack).get()) { trackIsSpokenFor=true; break; }
      }
      if (!trackIsSpokenFor) for (int iLep=0; iLep<gt->nLooseLep; iLep++) {
        if (!looseLeps[iLep]->matchedPF.isValid()) continue;
        if (softTrack==looseLeps[iLep]->matchedPF.get()) { trackIsSpokenFor=true; break; }
      }
      if (trackIsSpokenFor) continue;
      // Require tracks to have the lowest |dz| with the hardest PV amongst all others
      int idxVertexWithMinAbsDz=-1; float minAbsDz=9999;
      for (int iV=0; iV!=(int)event.vertices.size(); iV++) {
        auto& theVertex = event.vertices[iV];
        float vertexAbsDz = fabs(softTrack->dz(theVertex.position()));
        if (DEBUG > 10) 
          PDebug("PandaAnalyzer::JetHbbSoftActivity",Form("Track has |dz| %.2f with vertex %d",vertexAbsDz,iV));
        if (vertexAbsDz >= minAbsDz) continue;
        idxVertexWithMinAbsDz = iV;
        minAbsDz = vertexAbsDz;
      }
      if (idxVertexWithMinAbsDz!=0 || minAbsDz>0.2) continue;
      if (DEBUG > 10) 
        PDebug("PandaAnalyzer::JetHbbSoftActivity",
            Form("Track above 300 MeV has dz %.3f", 
              softTrack->track.isValid()?softTrack->track.get()->dz():-1));
      // Need to add High Quality track flags :-)
      bool trackIsInHbbEllipse=false; {
        double ellipse_x = softTrack->eta();
        double ellipse_y = softTrack->phi();
        double ellipse_term1 = pow(
          (TVector2::Phi_mpi_pi(ellipse_x - ellipse_h)*ellipse_cosA 
           + (ellipse_y - ellipse_k)*ellipse_sinA) / ellipse_a,
          2
        );
        double ellipse_term2 = pow(
          (TVector2::Phi_mpi_pi(ellipse_x - ellipse_h)*ellipse_sinA 
           - (ellipse_y - ellipse_k)*ellipse_cosA) / ellipse_b,
          2
        );
        double ellipse_equation = (ellipse_term1 + ellipse_term2);
        trackIsInHbbEllipse = (ellipse_equation <= 1.);
      } if (trackIsInHbbEllipse) continue;
      softTracksPJ.emplace_back(softTrack->px(),softTrack->py(),softTrack->pz(),softTrack->e());
    }
    if (DEBUG > 10) 
      PDebug("PandaAnalyzer::JetHbbSoftActivity",
          Form("Found %ld soft tracks that passed track quality cuts and the ellipse, jet constituency, and lepton matching vetoes",
            softTracksPJ.size()));
    softTrackJetDefinition = new fastjet::JetDefinition(fastjet::antikt_algorithm,0.4);
    fastjet::ClusterSequenceArea softTrackSequence(softTracksPJ, *softTrackJetDefinition, *areaDef);
    
    std::vector<fastjet::PseudoJet> softTrackJets(softTrackSequence.inclusive_jets(1.));
    if (DEBUG > 10) 
      PDebug("PandaAnalyzer::JetHbbSoftActivity",
          Form("Clustered %ld jets of pT>1GeV using anti-kT algorithm (dR 0.4) from the soft tracks",
            softTrackJets.size()));
    for (std::vector<fastjet::PseudoJet>::size_type iSTJ=0; iSTJ<softTrackJets.size(); iSTJ++) {
      if (fabs(softTrackJets[iSTJ].eta()) > 4.7) continue;
      gt->sumEtSoft1 += softTrackJets[iSTJ].Et(); 
      if (DEBUG > 10) 
        PDebug("PandaAnalyzer::JetHbbSoftActivity",
            Form("Soft jet %d has pT %.2f",
              (int)iSTJ,softTrackJets[iSTJ].pt()));
      if (softTrackJets[iSTJ].pt() >  2.)  gt->nSoft2++; else continue;
      if (softTrackJets[iSTJ].pt() >  5.)  gt->nSoft5++; else continue;
      if (softTrackJets[iSTJ].pt() > 10.) gt->nSoft10++; else continue;
    }
    tr->TriggerEvent("Soft activity");
  }

}
