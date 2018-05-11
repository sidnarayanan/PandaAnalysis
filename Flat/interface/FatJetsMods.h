#ifndef FATJETSMODS
#define FATJETSMODS

#include "Module.h"
#include "JetsMods.h" // need BaseJetMod 

namespace pa {
  class FatJetReclusterMod : public AnalysisMod {
  public: 
    FatJetReclusterMod(panda::EventAnalysis& event_, 
                       Config& cfg_,                 
                       Utils& utils_,                
                       GeneralTree& gt_,
                       int level_=0) :                 
      AnalysisMod("fjrecluster", event_, cfg_, utils_, gt_, level_) { 
        if (!on())
          return;
        jetDef.reset(new fastjet::JetDefinition(analysis.ak8 ? fastjet::antikt_algorithm : 
                                                               fastjet::cambridge_algorithm,
                                                analysis.ak8 ?  0.8 : 1.5));
    }
    virtual ~FatJetReclusterMod () { }

    virtual bool on() { return analysis.recluster; }
    
  protected:
    void do_init(Registry& registry) {
      fj1 = registry.accessConst<const panda::FatJet*>("fj1");
    }
    void do_execute();  
  private:
    std::shared_ptr<const panda::FatJet *const> fj1{nullptr}; // shared ptr to a const bare address
    std::unique_ptr<fastjet::JetDefinition> jetDef{nullptr};
  };


  class FatJetMod : public BaseJetMod {
  public: 
    FatJetMod(panda::EventAnalysis& event_, 
              Config& cfg_,                 
              Utils& utils_,                
              GeneralTree& gt_,
              int level_=0) :                 
      BaseJetMod("fatjet", event_, cfg_, utils_, gt_, level_),
      fj1(std::make_shared<const panda::FatJet*>(nullptr)),
      fatjets(analysis.ak8 ? event.puppiAK8Jets : event.puppiCA15Jets) { 
        recluster = addSubMod<FatJetReclusterMod>(); 
        jetType = "AK8PFPuppi";
    }
    virtual ~FatJetMod () { }

    virtual bool on() { return !analysis.genOnly && analysis.fatjet; }
    
  protected:
    void do_init(Registry& registry) {
      registry.publishConst("fj1", fj1);
      matchLeps = registry.accessConst<std::vector<panda::Lepton*>>("matchLeps");
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps"); 
      matchPhos = registry.accessConst<std::vector<panda::Photon*>>("tightPhos");
      dilep = registry.accessConst<TLorentzVector>("dilep"); 
    }
    void do_execute();  
    float getMSDCorr(float,float);
  private:
    void setupJES(); 

    std::shared_ptr<const panda::FatJet*> fj1{nullptr}; 
    panda::FatJetCollection &fatjets;

    std::shared_ptr<const std::vector<panda::Lepton*>> matchLeps{nullptr};
    std::shared_ptr<const std::vector<panda::Lepton*>> looseLeps{nullptr};
    std::shared_ptr<const std::vector<panda::Photon*>> matchPhos{nullptr};
    std::shared_ptr<const TLorentzVector> dilep{nullptr};

    FatJetReclusterMod *recluster{nullptr};
  };

  class FatJetMatchingMod : public AnalysisMod {
  public: 
    FatJetMatchingMod(panda::EventAnalysis& event_, 
                      Config& cfg_,                 
                      Utils& utils_,                
                      GeneralTree& gt_,
                      int level_=0) :                 
      AnalysisMod("fjmatching", event_, cfg_, utils_, gt_, level_) { }
    virtual ~FatJetMatchingMod () { }

    virtual bool on() { return !analysis.genOnly && analysis.fatjet && !analysis.isData; }
    
  protected:
    void do_init(Registry& registry) {
      fjPtr = registry.accessConst<const panda::FatJet*>("fj1");
      if (!analysis.isData)
        genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
    void do_execute();  
    void do_reset() { genObjects.clear(); }
  private:
    std::shared_ptr<const panda::FatJet *const> fjPtr{nullptr}; 
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr}; 

    std::map<const panda::GenParticle*,float> genObjects; // gen particle -> pt 
    const panda::GenParticle* matchGen(double eta, double phi, double r2, int pdgid=0) const;  
  };
}

#endif
