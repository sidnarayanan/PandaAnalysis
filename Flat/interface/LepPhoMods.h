#ifndef LEPPHOMODS
#define LEPPHOMODS

#include "Module.h"
#include "AnalyzerUtilities.h"
#include "array"
#include "PandaAnalysis/Utilities/interface/RoccoR.h"
#include "PandaAnalysis/Utilities/interface/RelIso.h"

namespace pa {

  class InclusiveLeptonMod : public AnalysisMod {
  public:
    InclusiveLeptonMod(panda::EventAnalysis& event_,
                       Config& cfg_,
                       Utils& utils_,
                       GeneralTree& gt_,
                       int level_=0) :
      AnalysisMod("inclep", event_, cfg_, utils_, gt_, level_),
      inclusiveLeps(std::make_shared<std::vector<panda::Lepton*>>()) { }
    virtual ~InclusiveLeptonMod () { }

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

  class SimpleLeptonMod : public AnalysisMod {
  public:
    SimpleLeptonMod(panda::EventAnalysis& event_,
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) :
      AnalysisMod("simplep", event_, cfg_, utils_, gt_, level_),
      looseLeps(std::v_make_shared<panda::Lepton*>()),
      matchLeps(std::v_make_shared<panda::Lepton*>()),
      lepPdgId(std::make_shared<std::array<int,4>>()),
      dilep(std::make_shared<TLorentzVector>()) { }
    virtual ~SimpleLeptonMod () { }

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

  class ComplicatedLeptonMod : public SimpleLeptonMod {
  public:
    ComplicatedLeptonMod(panda::EventAnalysis& event_,
                         Config& cfg_,
                         Utils& utils_,
                         GeneralTree& gt_,
                         int level_=0) :
      SimpleLeptonMod(event_, cfg_, utils_, gt_, level_) { name = "complep"; }
    virtual ~ComplicatedLeptonMod () { }

    virtual bool on() { return !analysis.genOnly && analysis.complicatedLeptons; }

  protected:
    void do_readData(TString dirPath);
    void do_execute();
    void do_init(Registry& registry) {
      SimpleLeptonMod::do_init(registry);
      if (!analysis.isData)
        genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
      if (analysis.hbb) 
        pfCandsMap = registry.access<EtaPhiMap<panda::PFCand>>("pfCandsMap"); 
    }
  private:
    std::unique_ptr<RoccoR> rochesterCorrection{nullptr};
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};
    std::shared_ptr<EtaPhiMap<panda::PFCand>> pfCandsMap{nullptr}; 
  };


  class SimplePhotonMod : public AnalysisMod {
  public:
    SimplePhotonMod(panda::EventAnalysis& event_,
                    Config& cfg_,
                    Utils& utils_,
                    GeneralTree& gt_,
                    int level_=0) :
      AnalysisMod("simplepho", event_, cfg_, utils_, gt_, level_),
      loosePhos(std::v_make_shared<panda::Photon*>()),
      tightPhos(std::v_make_shared<panda::Photon*>()) { }
    virtual ~SimplePhotonMod () { }

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
    void scaleFactors();
    std::shared_ptr<std::vector<panda::Photon*>> loosePhos, tightPhos;
  };


  class ComplicatedPhotonMod : public SimplePhotonMod {
  public:
    ComplicatedPhotonMod(panda::EventAnalysis& event_,
                         Config& cfg_,
                         Utils& utils_,
                         GeneralTree& gt_,
                         int level_=0) :
      SimplePhotonMod(event_, cfg_, utils_, gt_, level_) { name="comppho"; }
    virtual ~ComplicatedPhotonMod () { }

    virtual bool on() { return !analysis.genOnly && analysis.complicatedPhotons; }

  protected:
    void do_init(Registry& registry) {
      SimplePhotonMod::do_init(registry);
      matchLeps = registry.accessConst<std::vector<panda::Lepton*>>("matchLeps");
    }
    void do_execute();
  private:
    std::shared_ptr<const std::vector<panda::Lepton*>> matchLeps{nullptr};
    bool pfChargedPhotonMatch(const panda::Photon& photon);
  };


  class TauMod : public AnalysisMod {
  public:
    TauMod(panda::EventAnalysis& event_,
           Config& cfg_,
           Utils& utils_,
           GeneralTree& gt_,
           int level_=0) :
      AnalysisMod("tau", event_, cfg_, utils_, gt_, level_) { }
    virtual ~TauMod () { }

    virtual bool on() { return !analysis.genOnly; }

  protected:
    void do_init(Registry& registry) {
      matchLeps = registry.accessConst<std::vector<panda::Lepton*>>("matchLeps");
    }
    void do_execute();
  private:
    std::shared_ptr<const std::vector<panda::Lepton*>> matchLeps{nullptr};
  };

  class GenLepMod : public AnalysisMod {
  public:
    GenLepMod(panda::EventAnalysis& event_,
              Config& cfg_,
              Utils& utils_,
              GeneralTree& gt_,
              int level_=0) :
      AnalysisMod("genlep", event_, cfg_, utils_, gt_, level_) { }
    virtual ~GenLepMod () { }

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
