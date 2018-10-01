#ifndef HBBMODS
#define HBBMODS

#include "Module.h"
#include "PandaAnalysis/Utilities/interface/KinematicFit.h"

namespace pa {
  class GenVHMod : public AnalysisMod {
  public: 
    GenVHMod(panda::EventAnalysis& event_, 
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) : 
      AnalysisMod("genvh", event_, cfg_, utils_, gt_, level_) { }
    virtual ~GenVHMod () { }

    virtual bool on() { return !analysis.isData && analysis.vqqhbb; }
    
  protected:
    void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
    void do_execute();
  private:
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};
    template<typename LAMBDA, typename LAMBDB>
    void fillParticle(LAMBDA&& fn, LAMBDB&& childFn, float& pt, float& eta, float& phi, float& size) { 
      for (auto* pptr : *genP) {
        auto& p = pToGRef(pptr);
        if (!fn(abs(p.pdgid)))
          continue; 
        if (hasChild(p, *genP))
          continue; 
        pt = p.pt();
        eta = p.eta();
        phi = p.phi();
        size = -1; 
        // find children 
        for (auto* childptr : *genP) {
          auto& child = pToGRef(childptr); 
          if (!child.parent.isValid() || 
              child.parent.get() != pptr)
            continue; 
          if (!childFn(abs(child.pdgid)))
            continue; 
          size = std::max(static_cast<float>(DeltaR2(eta, phi, child.eta(), child.phi())), size);
        }
        return; 
      }
    }
  };

  class HbbMiscMod : public AnalysisMod {
  public: 
    HbbMiscMod(panda::EventAnalysis& event_, 
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) : 
      AnalysisMod("hbbmisc", event_, cfg_, utils_, gt_, level_) { }
    virtual ~HbbMiscMod () { }

    virtual bool on() { return !analysis.genOnly && analysis.hbb; }
    
  protected:
    void do_execute(); 
  };

  class KinFitMod : public AnalysisMod {
  public: 
    KinFitMod(panda::EventAnalysis& event_, 
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) : 
      AnalysisMod("zllhbbfit", event_, cfg_, utils_, gt_, level_) { 
      if (on()) {
        fit.reset(new kinfit::Fit(4,91));
        fit->setPrintLevel(-1); 
      }
    }
    virtual ~KinFitMod () { }

    virtual bool on() { return !analysis.genOnly && analysis.zllhbb; }
    
  protected:
    void do_init(Registry& registry) {
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
    }
    void do_execute(); 
  private:
    std::unique_ptr<kinfit::Fit> fit{nullptr}; 
    std::shared_ptr<const std::vector<panda::Lepton*>> looseLeps{nullptr};
  };

  class SoftActivityMod : public AnalysisMod {
  public: 
    SoftActivityMod(panda::EventAnalysis& event_, 
                    Config& cfg_,                 
                    Utils& utils_,                
                    GeneralTree& gt_,
                    int level_=0) :                 
      AnalysisMod("softactivity", event_, cfg_, utils_, gt_, level_) { }
    virtual ~SoftActivityMod () { }

    virtual bool on() { return !analysis.genOnly && analysis.hbb && !analysis.vbf && !analysis.vqqhbb && !analysis.bjetRegTraining; }
    
  protected:
    void do_init(Registry& registry) {
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
      jetDefSoftTrack.reset(new fastjet::JetDefinition(fastjet::antikt_algorithm,0.4));
    }
    void do_execute();  

  private:
    std::shared_ptr<std::vector<JESHandler>> jesShifts{nullptr}; 
    std::shared_ptr<const std::vector<panda::Lepton*>> looseLeps{nullptr};
    std::unique_ptr<fastjet::JetDefinition> jetDefSoftTrack{nullptr};
  };


  class GenJetNuMod : public AnalysisMod {
  public:
    GenJetNuMod(panda::EventAnalysis& event_, 
                 Config& cfg_,                 
                 Utils& utils_,                
                 GeneralTree& gt_,
                 int level_=0) :                 
      AnalysisMod("genjetnu", event_, cfg_, utils_, gt_, level_) { }
    virtual ~GenJetNuMod () { }

    bool on() { return !analysis.isData && analysis.bjetRegTraining && analysis.hbb; }
  protected:
    void do_init(Registry& registry) {
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
      jetDef.reset(new fastjet::JetDefinition(fastjet::antikt_algorithm,0.4));
    }
    void do_execute();
  private:
    std::shared_ptr<std::vector<JESHandler>> jesShifts{nullptr}; 
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};
    std::unique_ptr<fastjet::JetDefinition> jetDef{nullptr}; 
  };
}


#endif
