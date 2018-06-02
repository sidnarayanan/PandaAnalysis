#include "../interface/RelIso.h"

// Requires the relisomap to be filled with the event's PFCandidates
// Only use rho for muon calculations
double pa::MiniRelIso(
  panda::Lepton& lep, 
  EtaPhiMap<panda::PFCand>* pfCandsMapPtr, 
  decltype(panda::Event::rho) rho
) {
  // Get the parameter set
  auto abseta = std::abs(lep.eta());
  const auto& params = rho ? muonMiniIsoPars : (
    abseta < 1.4442 ? eleMiniIsoParsEB : eleMiniIsoParsEE
  );

  double chiso = 0;
  double nhiso = 0;
  double phiso = 0;
  double puiso = 0;

  // Copy of pat::miniIsoDr
  double drcut = std::max(params.min_dr, std::min(params.max_dr, params.kt_scale/lep.pt()));

  for (const auto* cand : pfCandsMapPtr->GetParticles(lep.eta(), lep.phi(), drcut)) {
    auto id = std::abs(cand->pdgId());
    auto pt = cand->pt();
    auto dr = sqrt(DeltaR2(lep.eta(), lep.phi(), cand->eta(), cand->phi()));
    switch(id) {
    case 211:
      if (cand->vertex.idx() == 0) {
        if (dr > params.deadcone_ch)
          chiso += pt;
      }
      else {
        if (pt > params.pt_threshold && dr > params.deadcone_pu)
          puiso += pt;
      }
      break;
    case 130:
      if (pt > params.pt_threshold && dr > params.deadcone_nh)
        nhiso += pt;
      break;
    case 22:
      if (pt > params.pt_threshold && dr > params.deadcone_ph)
        phiso += pt;
      break;
    // No action on default
    }
  }

  if (rho) {
    using ea_pair = std::pair<decltype(abseta), decltype(drcut)>;
    const std::vector<ea_pair> ea_corrections {
      {0.800, 0.0735},
      {1.300, 0.0619},
      {2.000, 0.0465},
      {2.200, 0.0433},
      {2.500, 0.0577}
    };
    auto ea_iter = std::upper_bound(ea_corrections.begin(), ea_corrections.end(), abseta,
                                    [] (decltype(abseta) val, const ea_pair& elem) { return val < elem.first; } );
    auto ea = ea_iter == ea_corrections.end() ? 0.0 : ea_iter->second;

    auto correction = rho * ea * std::pow(drcut/0.3, 2);
    return (chiso + std::max(0.0, nhiso + phiso - correction))/lep.pt();
  }

  return (chiso + std::max(0.0, nhiso + phiso - puiso))/lep.pt();
}

