#include "../interface/HbbOps.h"
#include "TVector2.h"

using namespace pa;
using namespace std;
using namespace panda;
namespace fj = fastjet;


void GenVHOp::do_execute()
{
  fillParticle([](int x)->bool { return x == 25; }, 
               [](int x)->bool { return x == 5; }, 
               gt.genHPt, gt.genHEta, gt.genHPhi, gt.genHSize);
  fillParticle([](int x)->bool { return x == 23 || x == 24; }, 
               [](int x)->bool { return x < 6; }, 
               gt.genVPt, gt.genVEta, gt.genVPhi, gt.genVSize);
}

void HbbMiscOp::do_execute()
{
  TVector2 trkmet;
  const auto& pv = event.vertices.at(0).position(); 
  for (auto& c : event.pfCandidates) {
    if (c.q() == 0 || fabs(c.dz(pv)) > 0.1) 
      continue;
    TVector2 v; v.SetMagPhi(c.pt(), c.phi());
    trkmet += v;
  }
  gt.trkmetDZ = trkmet.Mod();
  gt.trkmetDZphi = trkmet.Phi();
}

void KinFitOp::do_execute()
{
  if (gt.hbbm[0] <= 0 || looseLeps->size() < 2) {
    gt.hbbm_fit = gt.hbbm_dreg[0];
    return;
  }

  TLorentzVector res; res.SetPtEtaPhiM(gt.pfmet[0], 0, gt.pfmetphi[0], 0);
  std::array<TLorentzVector, 2> hbb;
  for (int i = 0; i != 2; ++i) {
    int idx = gt.hbbjtidx[0][i];
    hbb[i].SetPtEtaPhiM(gt.jotPt[0][idx],
                        gt.jotEta[idx],
                        gt.jotPhi[idx],
                        gt.jotM[idx]);
    hbb[i] *= gt.jotDeepBReg[i];
    res -= hbb[i];
    res -= (*looseLeps)[i]->p4(); 
  } 


  fit->setParticle(0, (*looseLeps)[0]->p4(), 0.01);
  fit->setParticle(1, (*looseLeps)[1]->p4(), 0.01);
  fit->setParticle(2, hbb[0], gt.jotDeepBRegWidth[0]);
  fit->setParticle(3, hbb[1], gt.jotDeepBRegWidth[1]);
  // fit->setParticle(4, res, res.Pt()/8); // 8 GeV res at 0 MET

  fit->run();

  gt.hbbm_fit = (hbb[0]*fit->getScale(2) + hbb[1]*fit->getScale(3)).M();
   
}

void SoftActivityOp::do_execute() 
{
  // Soft activity
  // TODO: I'm seeing a lot of sqrts replaced by pow(,2)s here. I think this can be sped up -SN
  int shift = jes2i(shiftjes::kNominal);
  auto& jets = (*jesShifts)[shift];

  if (gt.hbbm[shift] > 0.) {
    gt.sumEtSoft1=0; gt.nSoft2=0; gt.nSoft5=0; gt.nSoft10=0;
    // Define the ellipse of particles to forget about
    // https://math.stackexchange.com/questions/426150/
    // what-is-the-general-equation-of-the-ellipse-that-is-not-in-the-origin-and-rotate 
    // ((x-h)cos(A) + (y-k)sin(A))^2 /a^2 + ((x-h)sin(A) - (y-k)cos(A))^2 /b^2 <=1
    double ellipse_cosA, ellipse_sinA, ellipse_h, ellipse_k, ellipse_a, ellipse_b; {
      double ellipse_alpha;
      float phi1=gt.jotPhi[gt.hbbjtidx[shift][0]], phi2=gt.jotPhi[gt.hbbjtidx[shift][1]];
      float eta1=gt.jotEta[gt.hbbjtidx[shift][0]], eta2=gt.jotEta[gt.hbbjtidx[shift][1]];
      double phi1MinusPhi2 = phi1-phi2;
      double eta1MinusEta2 = eta1-eta2;
      double phi1MinusPhi2MPP = TVector2::Phi_mpi_pi(phi1MinusPhi2);
      ellipse_alpha = atan2( phi1MinusPhi2MPP, eta1MinusEta2);
      // compute delta R using already computed qty's to save time
      ellipse_a = (sqrt(pow(eta1MinusEta2,2) 
                   + pow(phi1MinusPhi2MPP,2)) + 1.)/2.; // Major axis 2*a = dR(b,b)+1
      ellipse_b = 1./2.; // Minor axis 2*b = 1
      ellipse_h = (eta1+eta2)/2.;
      ellipse_k = TVector2::Phi_mpi_pi(phi2 + phi1MinusPhi2MPP/2.);
      ellipse_cosA = cos(ellipse_alpha);
      ellipse_sinA = sin(ellipse_alpha);
    }

    // Find out which PF constituents to not use
    const RefVector<PFCand> &jet1Tracks = jets.cleaned[gt.hbbjtidx[shift][0]]->base->constituents,
                            &jet2Tracks = jets.cleaned[gt.hbbjtidx[shift][1]]->base->constituents;

    // Get vector of pseudo jets for clustering
    panda::PFCandCollection &allTracks = event.pfCandidates;
    vector<fj::PseudoJet> softTracksPJ;
    softTracksPJ.reserve(allTracks.size());
    panda::PFCand *softTrack=0;
    for (auto &softTrackRef : allTracks) {
      softTrack = &softTrackRef;
      // Minimum track pT threshold (300 MeV default)
      if (softTrack->pt() < cfg.minSoftTrackPt) 
        continue;
      if (fabs(softTrack->eta()) > 4.7)
        continue;
      // High quality track flag
      if (!softTrack->track.isValid() || !softTrack->track.get()->highPurity) 
        continue;
      // Only consider tracks with dz < 0.2 w.r.t. the primary vertex
      if (fabs(softTrack->track.get()->dz()) > 0.2) 
        continue;
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
      if (!trackIsSpokenFor) for (int iLep=0; iLep<gt.nLooseLep; iLep++) {
        if (!(*looseLeps)[iLep]->matchedPF.isValid()) continue;
        if (softTrack==(*looseLeps)[iLep]->matchedPF.get()) { trackIsSpokenFor=true; break; }
      }
      if (trackIsSpokenFor) 
        continue;
      // Require tracks to have the lowest |dz| with the hardest PV amongst all others
      int idxVertexWithMinAbsDz=-1; float minAbsDz=9999;
      for (int iV=0; iV!=(int)event.vertices.size(); iV++) {
        auto& theVertex = event.vertices[iV];
        float vertexAbsDz = fabs(softTrack->dz(theVertex.position()));
        if (vertexAbsDz >= minAbsDz) 
          continue;
        idxVertexWithMinAbsDz = iV;
        minAbsDz = vertexAbsDz;
      }
      if (idxVertexWithMinAbsDz!=0 || minAbsDz>0.2) 
        continue;
      // Need to add High Quality track flags :-)
      bool trackIsInHbbEllipse=false; {
        double ellipse_x = softTrack->eta();
        double ellipse_y = softTrack->phi();
        double ellipse_term1 = pow(
          ((ellipse_x - ellipse_h)*ellipse_cosA 
           + TVector2::Phi_mpi_pi(ellipse_y - ellipse_k)*ellipse_sinA) / ellipse_a,
          2
        );
        double ellipse_term2 = pow(
          ((ellipse_x - ellipse_h)*ellipse_sinA 
           - TVector2::Phi_mpi_pi(ellipse_y - ellipse_k)*ellipse_cosA) / ellipse_b,
          2
        );
        double ellipse_equation = (ellipse_term1 + ellipse_term2);
        trackIsInHbbEllipse = (ellipse_equation <= 1.);
      } if (trackIsInHbbEllipse) continue;
      softTracksPJ.emplace_back(softTrack->px(),softTrack->py(),softTrack->pz(),softTrack->e());
    }
    fj::ClusterSequenceArea softTrackSequence(softTracksPJ, *jetDefSoftTrack, *(utils.areaDef));
    
    vector<fj::PseudoJet> softTrackJets(softTrackSequence.inclusive_jets(1.));
    for (vector<fj::PseudoJet>::size_type iSTJ=0; iSTJ<softTrackJets.size(); iSTJ++) {
      if (fabs(softTrackJets[iSTJ].eta()) > 4.7) continue;
      gt.sumEtSoft1 += softTrackJets[iSTJ].Et(); 
      if (softTrackJets[iSTJ].pt() >  2.)  gt.nSoft2++; else continue;
      if (softTrackJets[iSTJ].pt() >  5.)  gt.nSoft5++; else continue;
      if (softTrackJets[iSTJ].pt() > 10.) gt.nSoft10++; else continue;
    }
  }
}

void GenJetNuOp::do_execute()
{
  vector<fj::PseudoJet> finalStates;
  vector<const panda::GenParticle*> bcs;
  for (auto* pptr : *genP) {
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

  fj::ClusterSequenceArea seq(finalStates, *jetDef, *(utils.areaDef));
  vector<fj::PseudoJet> allJets(seq.inclusive_jets(0.01));

  map<fj::PseudoJet*,int> flavorMap;
  for (auto &pj : allJets) {
    int flavor = 0;
    for (auto *bc : bcs) {
      if (DeltaR2(pj.eta(), pj.phi(), bc->eta(), bc->phi()) < 0.09) {
        flavor = abs(bc->pdgid);
        break;
      }
    }
    flavorMap[&pj] = flavor;
  }

  auto& jets = (*jesShifts)[0];
  unsigned N = jets.cleaned.size();
  for (unsigned i = 0; i != N; ++i) {
    const panda::Jet& reco = jets.cleaned[i]->get_base();
    for (auto &pj : allJets) {
      if (DeltaR2(pj.eta(), pj.phi(), reco.eta(), reco.phi()) < 0.09) {
        gt.jotGenPt[i] = pj.pt();
        gt.jotGenEta[i] = pj.eta();
        gt.jotGenPhi[i] = pj.phi();
        gt.jotGenM[i] = pj.m();
        gt.jotFlav[i] = flavorMap[&pj];
        gt.jotGenDEta[i] = gt.jotGenEta[i] - gt.jotEta[i];
        gt.jotGenDPhi[i] = SignedDeltaPhi(gt.jotGenPhi[i], gt.jotPhi[i]);
        break;
      }
    }
  }
}
