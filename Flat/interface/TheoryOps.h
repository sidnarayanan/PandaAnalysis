#ifndef THEORYOPS
#define THEORYOPS

#include "Operator.h"

namespace pa {
  class GenStudyEWKOp : public AnalysisOp {
  public: 
    GenStudyEWKOp(panda::EventAnalysis& event_, 
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) : 
      AnalysisOp("genstudyewk", event_, cfg_, utils_, gt_, level_) { }
    virtual ~GenStudyEWKOp () { }

    virtual bool on() { return !analysis.isData && !analysis.genOnly && (analysis.complicatedLeptons || analysis.complicatedPhotons); }
    
  protected:
    virtual void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
      loosePhos = registry.accessConst<std::vector<panda::Photon*>>("loosePhos");
      lepPdgId = registry.accessConst<std::array<int,4>>("lepPdgId"); 
    }
    virtual void do_execute(); 
  private:
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};
    std::shared_ptr<const std::vector<panda::Lepton*>> looseLeps{nullptr}; 
    std::shared_ptr<const std::vector<panda::Photon*>> loosePhos{nullptr}; 
    std::shared_ptr<const std::array<int,4>> lepPdgId{nullptr};
    
  };


  class QCDUncOp : public AnalysisOp {
  public: 
    QCDUncOp(panda::EventAnalysis& event_, 
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) : 
      AnalysisOp("qcdunc", event_, cfg_, utils_, gt_, level_) { }
    virtual ~QCDUncOp () { }

    virtual bool on() { return !analysis.isData; }
    
  protected:
    void do_execute(); 
  };


  class SignalGenOp : public AnalysisOp {
  public: 
    SignalGenOp(panda::EventAnalysis& event_, 
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) : 
      AnalysisOp("signalweight", event_, cfg_, utils_, gt_, level_) { }
    virtual ~SignalGenOp () { }

    virtual bool on() { return !analysis.isData && analysis.processType==kSignal; }
    
  protected:
    void do_execute(); 
    void do_init(Registry& registry) {
      wIDs = registry.accessConst<std::vector<TString>>("wIDs"); 
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
  private:
    std::shared_ptr<const std::vector<TString>> wIDs {nullptr};
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};
  };


  class HFCountingOp : public AnalysisOp {
  public: 
    HFCountingOp(panda::EventAnalysis& event_, 
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) : 
      AnalysisOp("hfcounting", event_, cfg_, utils_, gt_, level_) { }
    virtual ~HFCountingOp () { }

    virtual bool on() { return !analysis.isData; }
    
  protected:
    void do_execute(); 
    void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
  private:
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};
  };

  class KFactorOp : public AnalysisOp {
  public: 
    KFactorOp(panda::EventAnalysis& event_, 
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) : 
      AnalysisOp("kfactor", event_, cfg_, utils_, gt_, level_) { }
    virtual ~KFactorOp () { }

    virtual bool on() { return !analysis.isData; }
    
  protected:
    void do_execute(); 
    void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps", true);
    }
  private:
    void do_toppt(); 
    void do_vpt(); 

    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};
    std::shared_ptr<const std::vector<panda::Lepton*>> looseLeps{nullptr}; 
  };
}

#endif
