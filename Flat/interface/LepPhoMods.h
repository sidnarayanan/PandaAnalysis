#ifndef LEPPHOMODS
#define LEPPHOMODS

#include "Module.h"
#include "array"

namespace pa {
  class SimpleLeptonMod : public AnalysisMod {
  public: 
    SimpleLeptonMod(const panda::EventAnalysis& event_, 
                    const Config& cfg_,
                    const Utils& utils_,
                    GeneralTree& gt_) : 
      AnalysisMod("simplep", event_, cfg_, utils_, gt_) { }
    virtual ~SimpleLeptonMod () { }
    
  protected:
    virtual void do_initialize(Registry& registry) {
      registry.publishConst("looseLeps", &looseLeps);
      registry.publishConst("tightLeps", &tightLeps);
      registry.publishConst("lepPdgId", &lepPdgId);
    }
    virtual void do_execute(); 
    virtual void do_reset() { 
      looseLeps.clear();
      tightLeps.clear();
      for (auto& id : lepPdgId)
        id = 0;
    }
    void scaleFactors();

    std::vector<panda::Lepton*> looseLeps, tightLeps;
    std::array<int,4> lepPdgId;
    
  };


  class ComplicatedLeptonMod : public SimpleLeptonMod {
  public: 
    ComplicatedLeptonMod(const panda::EventAnalysis& event_, 
                         const Config& cfg_,
                         const Utils& utils_,
                         GeneralTree& gt_) : 
      AnalysisMod("complep", event_, cfg_, utils_, gt_) { }
    ~ComplicatedLeptonMod () { delete rochesterCorrection; }
    
  protected:
    void do_readData(TString dirPath);
    void do_execute(); 
  private:
    RoccoR *rochesterCorrection{nullptr};
  };

  class InclusiveLeptonMod : public AnalysisMod {
  public: 
    InclusiveLeptonMod(const panda::EventAnalysis& event_, 
                       const Config& cfg_,
                       const Utils& utils_,
                       GeneralTree& gt_) : 
      AnalysisMod("inclep", event_, cfg_, utils_, gt_) { }
    virtual ~InclusiveLeptonMod () { }
    
  protected:
    virtual void do_initialize(Registry& registry) {
      registry.publishConst("inclusiveLeps", &inclusiveLeps);
    }
    virtual void do_execute(); 
    virtual void do_reset() { 
      inclusiveLeps.clear();
    }
  private:
    std::vector<panda::Lepton*> inclusiveLeps;
  };


  class SimplePhotonMod : public AnalysisMod {
  public: 
    SimplePhotonMod(const panda::EventAnalysis& event_, 
                    const Config& cfg_,
                    const Utils& utils_,
                    GeneralTree& gt_) : 
      AnalysisMod("simplepho", event_, cfg_, utils_, gt_) { }
    virtual ~SimplePhotonMod () { }
    
  protected:
    virtual void do_initialize(Registry& registry) {
      registry.publishConst("loosePhos", &loosePhos);
    }
    virtual void do_execute(); 
    virtual void do_reset() { 
      loosePhos.clear();
    }
    void scaleFactors();
    std::vector<panda::Photon*> loosePhos;
  };


  class ComplicatedPhotonMod : public SimplePhotonMod {
  public: 
    ComplicatedPhotonMod(const panda::EventAnalysis& event_, 
                         const Config& cfg_,
                         const Utils& utils_,
                         GeneralTree& gt_) : 
      AnalysisMod("comppho", event_, cfg_, utils_, gt_) { }
    ~ComplicatedPhotonMod () { }
    
  protected:
    void do_initalize(Registry& registry) {
      SimplePhotonMod::do_initialize(registry);
      matchLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
    }
    void do_execute(); 
  private:
    const std::vector<panda::Lepton*> *matchLeps; 
    bool pfChargedPhotonMatch(const panda::Photon& photon);
  };


  class TauMod : public AnalysisMod {
  public: 
    TauMod(const panda::EventAnalysis& event_, 
           const Config& cfg_,
           const Utils& utils_,
           GeneralTree& gt_) : 
      AnalysisMod("tau", event_, cfg_, utils_, gt_) { }
    ~TauMod () { }
    
  protected:
    void do_initalize(Registry& registry) {
      matchLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
    }
    void do_execute(); 
  private:
    const std::vector<panda::Lepton*> *matchLeps; 
  };

  class GenLepMod : public AnalysisMod {
  public: 
    GenLepMod(const panda::EventAnalysis& event_, 
              const Config& cfg_,
              const Utils& utils_,
              GeneralTree& gt_) : 
      AnalysisMod("genlep", event_, cfg_, utils_, gt_) { }
    ~GenLepMod () { }
    
  protected:
    void do_execute(); 
  };

}

#endif
