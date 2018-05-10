#ifndef COMMONMODS
#define COMMONMODS

#include "Module.h"
#include "AnalyzerUtilities.h"

namespace pa {
  class RecoilMod : public AnalysisMod {
  public:
    RecoilMod(panda::EventAnalysis& event_, 
              Config& cfg_,                 
              Utils& utils_,                
              GeneralTree& gt_,
              int level_=0) :                 
      AnalysisMod("recoil", event_, cfg_, utils_, gt_, level_) { }
    virtual ~RecoilMod () { }

    bool on() { return !analysis.genOnly && analysis.recoil; }

  protected:
    void do_init(Registry& registry) {
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
      loosePhos = registry.accessConst<std::vector<panda::Photon*>>("loosePhos");
      lepPdgId = registry.accessConst<std::array<int,4>>("lepPdgId");
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
    }
    void do_execute();

  private:
    const std::vector<panda::Lepton*> *looseLeps{nullptr};
    const std::vector<panda::Photon*> *loosePhos{nullptr};
    const std::array<int,4> *lepPdgId {nullptr};
    std::vector<JESHandler> *jesShifts{nullptr};
  };

  class TriggerMod : public AnalysisMod {
  public:
    TriggerMod(panda::EventAnalysis& event_, 
               Config& cfg_,                 
               Utils& utils_,                
               GeneralTree& gt_,
               int level_=0) :                 
      AnalysisMod("trigger", event_, cfg_, utils_, gt_, level_),
      triggerHandlers(kNTrig) { }
    virtual ~TriggerMod () { }

    bool on() { return !analysis.genOnly && (analysis.isData || analysis.mcTriggers); }

  protected:
    void do_init(Registry& registry);
    void do_execute();

  private:
    void checkEle32();
    std::vector<TriggerHandler> triggerHandlers; 
  };

  class TriggerEffMod : public AnalysisMod {
  public:
    TriggerEffMod(panda::EventAnalysis& event_, 
               Config& cfg_,                 
               Utils& utils_,                
               GeneralTree& gt_,
               int level_=0) :                 
      AnalysisMod("triggereff", event_, cfg_, utils_, gt_, level_) { }
    virtual ~TriggerEffMod () { }

    bool on() { return !analysis.genOnly && !analysis.isData; }

  protected:
    void do_init(Registry& registry) {
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps"); 
    }
    void do_execute();

  private:
    const std::vector<panda::Lepton*> *looseLeps{nullptr};
  };


  class GlobalMod : public AnalysisMod {
  public:
    GlobalMod(panda::EventAnalysis& event_, 
               Config& cfg_,                 
               Utils& utils_,                
               GeneralTree& gt_,
               int level_=0) :                 
      AnalysisMod("global", event_, cfg_, utils_, gt_, level_),
      jesShifts(jes2i(shiftjes::N)) { 
        JESLOOP {
          jesShifts[shift].shift_idx = shift;
        }
      }
    virtual ~GlobalMod () { }
    
  protected:
    void do_init(Registry& registry) { registry.publish("jesShifts", &jesShifts); }
    void do_execute();  
    void do_reset() { 
      for (auto& s : jesShifts)
        s.clear();
    }
  private:
    std::vector<JESHandler> jesShifts;
  };

  class GenPMod : public AnalysisMod {
  public:
    GenPMod(panda::EventAnalysis& event_, 
            Config& cfg_,                 
            Utils& utils_,                
            GeneralTree& gt_,
            int level_=0) :                 
      AnalysisMod("gendup", event_, cfg_, utils_, gt_, level_) { }
    virtual ~GenPMod () { }

    bool on() { return !analysis.isData; }
    
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
    void merge_particles(panda::Collection<T>& genParticles) {
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
