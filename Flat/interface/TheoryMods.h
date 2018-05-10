#ifndef THEORYMODS
#define THEORYMODS

#include "Module.h"

namespace pa {
  class GenStudyEWKMod : public AnalysisMod {
  public: 
    GenStudyEWKMod(panda::EventAnalysis& event_, 
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) : 
      AnalysisMod("genstudyewk", event_, cfg_, utils_, gt_, level_) { }
    virtual ~GenStudyEWKMod () { }

    virtual bool on() { return !analysis.isData && (analysis.complicatedLeptons || analysis.complicatedPhotons); }
    
  protected:
    virtual void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
      loosePhos = registry.accessConst<std::vector<panda::Photon*>>("loosePhos");
      lepPdgId = registry.accessConst<std::array<int,4>>("lepPdgId"); 
    }
    virtual void do_execute(); 
  private:
    const std::vector<panda::Particle*> *genP{nullptr};
    const std::vector<panda::Lepton*> *looseLeps{nullptr}; 
    const std::vector<panda::Photon*> *loosePhos{nullptr}; 
    const std::array<int,4> *lepPdgId{nullptr};
    
  };


  class QCDUncMod : public AnalysisMod {
  public: 
    QCDUncMod(panda::EventAnalysis& event_, 
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) : 
      AnalysisMod("qcdunc", event_, cfg_, utils_, gt_, level_) { }
    virtual ~QCDUncMod () { }

    virtual bool on() { return !analysis.isData; }
    
  protected:
    void do_execute(); 
  };


  class SignalGenMod : public AnalysisMod {
  public: 
    SignalGenMod(panda::EventAnalysis& event_, 
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) : 
      AnalysisMod("signalweight", event_, cfg_, utils_, gt_, level_) { }
    virtual ~SignalGenMod () { }

    virtual bool on() { return !analysis.isData && analysis.processType==kSignal; }
    
  protected:
    void do_execute(); 
    void do_init(Registry& registry) {
      wIDs = registry.accessConst<std::vector<TString>>("wIDs"); 
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
  private:
    const std::vector<TString> *wIDs {nullptr};
    const std::vector<panda::Particle*> *genP{nullptr};
  };


  class HFCountingMod : public AnalysisMod {
  public: 
    HFCountingMod(panda::EventAnalysis& event_, 
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) : 
      AnalysisMod("hfcounting", event_, cfg_, utils_, gt_, level_) { }
    virtual ~HFCountingMod () { }

    virtual bool on() { return !analysis.isData; }
    
  protected:
    void do_execute(); 
    void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
  private:
    const std::vector<panda::Particle*> *genP{nullptr};
  };

  class KFactorMod : public AnalysisMod {
  public: 
    KFactorMod(panda::EventAnalysis& event_, 
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) : 
      AnalysisMod("kfactor", event_, cfg_, utils_, gt_, level_) { }
    virtual ~KFactorMod () { }

    virtual bool on() { return !analysis.isData; }
    
  protected:
    void do_execute(); 
    void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
    }
  private:
    void toppt(); 
    void vpt(); 

    const std::vector<panda::Particle*> *genP{nullptr};
    const std::vector<panda::Lepton*> *looseLeps{nullptr}; 
  };
}

#endif
