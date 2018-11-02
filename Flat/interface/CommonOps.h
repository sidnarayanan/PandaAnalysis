#ifndef COMMONOPS
#define COMMONOPS

#include "Operator.h"
#include "AnalyzerUtilities.h"
#include "PandaAnalysis/Utilities/interface/EtaPhiMap.h"

namespace pa {
  class MapOp : public AnalysisOp {
  public:
    MapOp(panda::EventAnalysis& event_,
              Config& cfg_,
              Utils& utils_,
              GeneralTree& gt_,
              int level_=0) :
      AnalysisOp("map", event_, cfg_, utils_, gt_, level_) { }
    virtual ~MapOp () { }

    bool on() { return analysis.complicatedLeptons; }

  protected:
    void do_init(Registry& registry) {
      if (analysis.hbb) {
        pfCandsMap = std::make_shared<EtaPhiMap<panda::PFCand>>(0.1);
        registry.publish("pfCandsMap", pfCandsMap);
      }
      if (!analysis.isData && false) { // figure out if this is useful first
        genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
        genPMap = std::make_shared<EtaPhiMap<panda::Particle*>>(0.3, 5.0,
                                    [] (panda::Particle*const* p)->double { return (*p)->eta(); },
                                    [] (panda::Particle*const* p)->double { return (*p)->phi(); });
        registry.publish("genPMap", genPMap);
      }
    }
    void do_execute() {
      if (pfCandsMap.get() != nullptr) {
        pfCandsMap->AddParticles(event.pfCandidates);
      }
      if (genPMap.get() != nullptr) {
        genPMap->AddParticles(*(genP.get()));
      }
    }

  private:
    std::shared_ptr<EtaPhiMap<panda::PFCand>> pfCandsMap{nullptr}; 
    std::shared_ptr<EtaPhiMap<panda::Particle*>> genPMap{nullptr}; 
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};
  };

  class RecoilOp : public AnalysisOp {
  public:
    RecoilOp(panda::EventAnalysis& event_,
              Config& cfg_,
              Utils& utils_,
              GeneralTree& gt_,
              int level_=0) :
      AnalysisOp("recoil", event_, cfg_, utils_, gt_, level_) { }
    virtual ~RecoilOp () { }

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
    std::shared_ptr<const std::vector<panda::Lepton*>> looseLeps{nullptr};
    std::shared_ptr<const std::vector<panda::Photon*>> loosePhos{nullptr};
    std::shared_ptr<const std::array<int,4>> lepPdgId {nullptr};
    std::shared_ptr<std::vector<JESHandler>> jesShifts{nullptr};
  };

  class TriggerOp : public AnalysisOp {
  public:
    TriggerOp(panda::EventAnalysis& event_,
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) :
      AnalysisOp("trigger", event_, cfg_, utils_, gt_, level_),
      triggerHandlers(kNTrig) { }
    virtual ~TriggerOp () { }

    bool on() { return !analysis.genOnly && (analysis.isData || analysis.mcTriggers); }

  protected:
    void do_init(Registry& registry);
    void do_execute();

  private:
    void checkEle32();
    std::vector<TriggerHandler> triggerHandlers;
  };

  class TriggerEffOp : public AnalysisOp {
  public:
    TriggerEffOp(panda::EventAnalysis& event_,
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) :
      AnalysisOp("triggereff", event_, cfg_, utils_, gt_, level_) { }
    virtual ~TriggerEffOp () { }

    bool on() { return !analysis.genOnly && !analysis.isData; }

  protected:
    void do_init(Registry& registry) {
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
    }
    void do_execute();

  private:
    std::shared_ptr<const std::vector<panda::Lepton*>> looseLeps{nullptr};
  };


  class GlobalOp : public AnalysisOp {
  public:
    GlobalOp(panda::EventAnalysis& event_,
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) :
      AnalysisOp("global", event_, cfg_, utils_, gt_, level_),
      jesShifts(std::v_make_shared<JESHandler>(jes2i(shiftjes::N))) {
        JESLOOP {
          (*jesShifts)[shift].shift_idx = shift;
        }
      }
    virtual ~GlobalOp () { }

  protected:
    void do_init(Registry& registry) {
      registry.publish("jesShifts", jesShifts);
      auto dummy = registry.access<std::vector<JESHandler>>("jesShifts");
    }
    void do_execute();
    void do_reset() {
      for (auto& s : *jesShifts)
        s.clear();
    }
  private:
    std::shared_ptr<std::vector<JESHandler>> jesShifts;
  };

  template <typename T>
  class BaseGenPOp : public BaseAnalysisOp<T> {
  public:
    BaseGenPOp(panda::EventAnalysis& event_,
                Config& cfg_,
                Utils& utils_,
                T& gt_,
                int level_=0) :
      BaseAnalysisOp<T>("gendup", event_, cfg_, utils_, gt_, level_),
      genP(std::v_make_shared<panda::Particle*>()) { }
    virtual ~BaseGenPOp () { }

    bool on() { return !this->analysis.isData; }

  protected:
    void do_init(Registry& registry) {
      registry.publishConst("genP", genP);
    }
    void do_execute();
    void do_reset() {
      genP->clear();
    }
  private:
    std::shared_ptr<std::vector<panda::Particle*>> genP;

    template <typename P>
    void merge_particles(panda::Collection<P>& genParticles) {
      genP->reserve(genParticles.size()); // approx
      for (auto& g : genParticles) {
        bool foundDup = false;
        if (g.finalState) {
          float ptThreshold = g.pt() * 0.01;
          for (auto* pPtr : *genP) {
            const P* gPtr = static_cast<const P*>(pPtr);
            if (!gPtr->finalState)
              continue;
            if ((g.pdgid == gPtr->pdgid) &&
                (fabs(g.pt() - gPtr->pt()) < ptThreshold) &&
                (DeltaR2(g.eta(), g.phi(), gPtr->eta(), gPtr->phi()) < 0.001)) {
              foundDup = true;
              if (this->cfg.DEBUG > 8) {
                logger.debug("Found duplicate",
                       Form("p1(%8.3f,%5.1f,%5.1f,%5i,%i) <-> p2(%8.3f,%5.1f,%5.1f,%5i,%i)",
                            g.pt(), g.eta(), g.phi(), g.pdgid, g.finalState ? 1 : 0,
                            gPtr->pt(), gPtr->eta(), gPtr->phi(), gPtr->pdgid, gPtr->finalState ? 1 : 0));
              }
              break;
            } // matches
          } // genP loop
        } // if final state
        if (!foundDup) {
          genP->push_back(&g);
        }
      }
    }
  };
  typedef BaseGenPOp<GeneralTree> GenPOp;
  typedef BaseGenPOp<HeavyResTree> HRGenPOp;
}

#endif
