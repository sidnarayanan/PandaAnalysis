#ifndef DEEPMODS
#define DEEPMODS

#include "Module.h"
#include "set"
#include "PandaAnalysis/Utilities/interface/EnergyCorrelations.h"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Njettiness.hh"

namespace pa {
  template <typename GENP>
  class DeepGenMod : public AnalysisMod {
  public: 
    DeepGenMod(const panda::EventAnalysis& event_, 
               const Config& cfg_,
               const Utils& utils_,
               GeneralTree& gt_);
    ~DeepGenMod () { 
      delete jetDef; delete softDrop; delete tauN;
      delete ecfcalc; delete grid;
    }
    
  protected:
    void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
      incrementAux(false); 
    }
    void do_execute(); 
    void do_reset() { 
      genJetInfo.reset();
      if (grid != nullptr)
        grid->clear();
    }
    void do_terminate() {
      incrementAux(true);
    }
    void countGenPartons(std::unordered_set<const GENP*>&);
    void incrementAux(bool close = false);
  private:
    const std::vector<panda::Particle*> *genP{nullptr};
    fastjet::JetDefinition       *jetDef                {nullptr};
    fastjet::contrib::SoftDrop   *softDrop              {nullptr};
    fastjet::contrib::Njettiness *tauN                  {nullptr};
    ECFCalculator *ecfcalc{nullptr};
    ParticleGridder *grid{nullptr};
    GenJetInfo genJetInfo;
    TFile *fAux{nullptr};
    TTree *tAux{nullptr};
    int auxCounter{0}; 
  };
}

#endif
