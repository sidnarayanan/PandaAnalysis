#ifndef LEPPHOOPS
#define LEPPHOOPS

#include "Operator.h"
#include "AnalyzerUtilities.h"
#include "array"
#include "PandaAnalysis/Utilities/interface/RoccoR.h"
#include "PandaAnalysis/Utilities/interface/RelIso.h"

namespace pa {

  class InclusiveLeptonOp : public AnalysisOp {
  public:
    InclusiveLeptonOp(panda::EventAnalysis& event_,
                       Config& cfg_,
                       Utils& utils_,
                       GeneralTree& gt_,
                       int level_=0) :
      AnalysisOp("inclep", event_, cfg_, utils_, gt_, level_),
      inclusiveLeps(std::make_shared<std::vector<panda::Lepton*>>()) { }
    virtual ~InclusiveLeptonOp () { }

    virtual bool on() { return !analysis.genOnly && analysis.hbb; }

  protected:
    virtual void do_init(Registry& registry) {
      registry.publishConst("inclusiveLeps", inclusiveLeps);
    }
    virtual void do_execute();
    virtual void do_reset() {
      inclusiveLeps->clear();
    }
  private:
    std::shared_ptr<std::vector<panda::Lepton*>> inclusiveLeps;
  };

  class SimpleLeptonOp : public AnalysisOp {
  public:
    SimpleLeptonOp(panda::EventAnalysis& event_,
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) :
      AnalysisOp("simplep", event_, cfg_, utils_, gt_, level_),
      looseLeps(std::v_make_shared<panda::Lepton*>()),
      matchLeps(std::v_make_shared<panda::Lepton*>()),
      lepPdgId(std::make_shared<std::array<int,4>>()),
      dilep(std::make_shared<TLorentzVector>()) { }
    virtual ~SimpleLeptonOp () { }

    virtual bool on() { return !analysis.genOnly && !analysis.complicatedLeptons; }

  protected:
    virtual void do_init(Registry& registry) {
      registry.publishConst("looseLeps", looseLeps); // sink sheps
      registry.publishConst("matchLeps", matchLeps);
      registry.publishConst("lepPdgId", lepPdgId);
      registry.publishConst("dilep", dilep);
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
    }
    virtual void do_execute();
    virtual void do_reset() {
      looseLeps->clear();
      matchLeps->clear();
      for (auto& id : *lepPdgId)
        id = 0;
      dilep->SetPtEtaPhiM(0,0,0,0);
    }
    void scaleFactors();

    std::shared_ptr<std::vector<panda::Lepton*>> looseLeps, matchLeps;
    std::shared_ptr<std::array<int,4>> lepPdgId;
    std::shared_ptr<std::vector<JESHandler>> jesShifts{nullptr};
    std::shared_ptr<TLorentzVector> dilep;

  };

  class ComplicatedLeptonOp : public SimpleLeptonOp {
  public:
    ComplicatedLeptonOp(panda::EventAnalysis& event_,
                         Config& cfg_,
                         Utils& utils_,
                         GeneralTree& gt_,
                         int level_=0) :
      SimpleLeptonOp(event_, cfg_, utils_, gt_, level_) { name = "complep"; }
    virtual ~ComplicatedLeptonOp () { }

    virtual bool on() { return !analysis.genOnly && analysis.complicatedLeptons; }

  protected:
    void do_readData(TString dirPath);
    void do_execute();
    void do_init(Registry& registry) {
      SimpleLeptonOp::do_init(registry);
      if (analysis.hbb) 
        pfCandsMap = registry.access<EtaPhiMap<panda::PFCand>>("pfCandsMap"); 
    }
  private:
    std::unique_ptr<RoccoR> rochesterCorrection{nullptr};
    std::shared_ptr<EtaPhiMap<panda::PFCand>> pfCandsMap{nullptr}; 
  };


  class SimplePhotonOp : public AnalysisOp {
  public:
    SimplePhotonOp(panda::EventAnalysis& event_,
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) :
      AnalysisOp("simplepho", event_, cfg_, utils_, gt_, level_),
      loosePhos(std::v_make_shared<panda::Photon*>()),
      tightPhos(std::v_make_shared<panda::Photon*>()) { }
    virtual ~SimplePhotonOp () { }

    virtual bool on() { return !analysis.genOnly && !analysis.complicatedPhotons; }

  protected:
    virtual void do_init(Registry& registry) {
      registry.publishConst("loosePhos", loosePhos);
      registry.publishConst("tightPhos", tightPhos);
    }
    virtual void do_execute();
    virtual void do_reset() {
      loosePhos->clear();
      tightPhos->clear();
    }
    bool is_tight(panda::Photon& p ) { // cannot use pointer to member reference
      return analysis.vbfhbb ?  p.tight : p.medium; 
    }
    void scaleFactors();
    std::shared_ptr<std::vector<panda::Photon*>> loosePhos, tightPhos;
  };


  class ComplicatedPhotonOp : public SimplePhotonOp {
  public:
    ComplicatedPhotonOp(panda::EventAnalysis& event_,
                         Config& cfg_,
                         Utils& utils_,
                         GeneralTree& gt_,
                         int level_=0) :
      SimplePhotonOp(event_, cfg_, utils_, gt_, level_) { name="comppho"; }
    virtual ~ComplicatedPhotonOp () { }

    virtual bool on() { return !analysis.genOnly && analysis.complicatedPhotons; }

  protected:
    void do_init(Registry& registry) {
      SimplePhotonOp::do_init(registry);
      matchLeps = registry.accessConst<std::vector<panda::Lepton*>>("matchLeps");
    }
    void do_execute();
  private:
    std::shared_ptr<const std::vector<panda::Lepton*>> matchLeps{nullptr};
    bool pfChargedPhotonMatch(const panda::Photon& photon);
  };


  class TauOp : public AnalysisOp {
  public:
    TauOp(panda::EventAnalysis& event_,
           Config& cfg_,
           Utils& utils_,
           GeneralTree& gt_,
           int level_=0) :
      AnalysisOp("tau", event_, cfg_, utils_, gt_, level_) { }
    virtual ~TauOp () { }

    virtual bool on() { return !analysis.genOnly; }

  protected:
    void do_init(Registry& registry) {
      matchLeps = registry.accessConst<std::vector<panda::Lepton*>>("matchLeps");
    }
    void do_execute();
  private:
    std::shared_ptr<const std::vector<panda::Lepton*>> matchLeps{nullptr};
  };

  class GenLepOp : public AnalysisOp {
  public:
    GenLepOp(panda::EventAnalysis& event_,
              Config& cfg_,
              Utils& utils_,
              GeneralTree& gt_,
              int level_=0) :
      AnalysisOp("genlep", event_, cfg_, utils_, gt_, level_) { }
    virtual ~GenLepOp () { }

    virtual bool on() { return analysis.vbf && !analysis.isData; }

  protected:
    void do_execute();
    void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
  private:
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};
  };

}

#endif
