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
