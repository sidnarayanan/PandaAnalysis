#ifndef RELISO 
#define RELISO 

#include "PandaAnalysis/Utilities/interface/EtaPhiMap.h"
#include "PandaTree/Objects/interface/Event.h"
#include "PandaCore/Tools/interface/Common.h"

namespace pa {

  class MiniIsoParams { // From Dan
  public:
    using type = double;
    type min_dr;
    type max_dr;
    type kt_scale;
    type pt_threshold;
    type deadcone_ch;
    type deadcone_pu;
    type deadcone_ph;
    type deadcone_nh;
  };
  const MiniIsoParams muonMiniIsoPars  {0.05, 0.2, 10.0, 0.5, 0.0001, 0.010, 0.01, 0.01};
  const MiniIsoParams eleMiniIsoParsEE {0.05, 0.2, 10.0, 0.0, 0.0150, 0.015, 0.08, 0.00};
  const MiniIsoParams eleMiniIsoParsEB {0.05, 0.2, 10.0, 0.0, 0.0000, 0.000, 0.00, 0.00};
  
  double MiniRelIso(
    panda::Lepton& lep, 
    EtaPhiMap<panda::PFCand>* pfCandsMapPtr, 
    decltype(panda::Event::rho) rho = 0
  );
}

#endif
