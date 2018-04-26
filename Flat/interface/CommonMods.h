#ifndef COMMONMODS
#define COMMONMODS

#include "Module.cc"
#include "AnalyzerUtilities.h"

namespace pa {
  class RecoilMod : public AnalysisMod {
  public:
    RecoilMod(const panda::EventAnalysis& event_, 
              const Config& cfg_,                 
              const Utils& utils_,                
              GeneralTree& gt_) :                 
      AnalysisMod("recoil", event_, cfg_, utils_, gt_) { }
    ~RecoilMod () { }

  protected:
    void do_initalize(Registry& registry) {
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
      loosePhos = registry.accessConst<std::vector<panda::Photon*>>("loosePhos");
      lepPdgId = registry.accessConst<std::array<int,4>>("lepPdgId");
    }
    void do_execute();

  private:
    const std::vector<panda::Lepton*> *looseLeps{nullptr};
    const std::vector<panda::Photon*> *loosePhos{nullptr};
    const std::array<int,4> lepPdgId *{nullptr};
  };

  class TriggerMod : public AnalysisMod {
  public:
    TriggerMod(const panda::EventAnalysis& event_, 
               const Config& cfg_,                 
               const Utils& utils_,                
               GeneralTree& gt_) :                 
      AnalysisMod("trigger", event_, cfg_, utils_, gt_) { }
    ~TriggerMod () { }

  protected:
    void do_initalize(Registry& registry);
    void do_execute();

  private:
    auto triggerHandlers = std::vector<TriggerHandler>(kNTrig) 
  };

  class GlobalMod : public AnalysisMod {
  public:
    GlobalMod(const panda::EventAnalysis& event_, 
               const Config& cfg_,                 
               const Utils& utils_,                
               GeneralTree& gt_) :                 
      AnalysisMod("global", event_, cfg_, utils_, gt_) { }
    ~GlobalMod () { }
    
  protected:
    void do_init(Registry& registry) { registry.publish("jesShifts", &jesShifts); }
    void do_execute();  
    void do_reset() { 
      for (auto& s : jesShifts)
        s.reset();
    }
  private:
    auto jesShifts = std::vector<JESHandler>(jes2i(shiftjes::N)); 
  };

  class GenPMod : public AnalysisMod {
  public:
    GenPMod(const panda::EventAnalysis& event_, 
            const Config& cfg_,                 
            const Utils& utils_,                
            GeneralTree& gt_) :                 
      AnalysisMod("gendup", event_, cfg_, utils_, gt_) { }
    ~GenPMod () { }
    
  protected:
    void do_init(Registry& registry) {
      registry.publishConst("genP", &genP);
    }
    void do_execute();  
    void do_reset() {
      genP.clear();
    }
  private:
    std::vector<panda::Particle*> genP;

    template <typename T>
    void merge_particles(const panda::Collection<T>& genParticles) {
      genP.reserve(genParticles.size()); // approx
      for (auto& g : genParticles) {
        bool foundDup = false;
        if (g.finalState) {
          float ptThreshold = g.pt() * 0.01;
          for (auto* pPtr : genP) {
            const T* gPtr = static_cast<const T*>(pPtr);
            if (!gPtr->finalState)
              continue;
            if ((g.pdgid == gPtr->pdgid) &&
                (fabs(g.pt() - gPtr->pt()) < ptThreshold) &&
                (DeltaR2(g.eta(), g.phi(), gPtr->eta(), gPtr->phi()) < 0.001)) {
              foundDup = true;
              if (cfg.DEBUG > 8) {
                PDebug("Found duplicate",
                       Form("p1(%8.3f,%5.1f,%5.1f,%5i,%i) <-> p2(%8.3f,%5.1f,%5.1f,%5i,%i)",
                            g.pt(), g.eta(), g.phi(), g.pdgid, g.finalState ? 1 : 0,
                            gPtr->pt(), gPtr->eta(), gPtr->phi(), gPtr->pdgid, gPtr->finalState ? 1 : 0));
              }
              break;
            } // matches
          } // genP loop
        } // if final state
        if (!foundDup) {
          genP.push_back(&g);
        }
      } 
    }
  };
}

#endif
