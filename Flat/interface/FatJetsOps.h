#ifndef FATJETSOPS
#define FATJETSOPS

#include "Operator.h"
#include "JetsOps.h" // need BaseJetOp
#include "PandaAnalysis/Utilities/interface/EnergyCorrelations.h"
#include "PandaAnalysis/Utilities/interface/HEPTopTaggerWrapperV2.h"
#include "fastjet/contrib/Njettiness.hh"

namespace pa {
  class FatJetReclusterOp : public AnalysisOp {
  public:
    FatJetReclusterOp(panda::EventAnalysis& event_,
                       Config& cfg_,
                       Utils& utils_,
                       GeneralTree& gt_,
                       int level_=0) :
      AnalysisOp("fjrecluster", event_, cfg_, utils_, gt_, level_) {
        if (!on())
          return;
        jetDef.reset(new fastjet::JetDefinition(analysis.ak8 ? fastjet::antikt_algorithm :
                                                               fastjet::cambridge_algorithm,
                                                analysis.ak8 ?  0.8 : 1.5));
    }
    virtual ~FatJetReclusterOp () { }

    virtual bool on() { return analysis.recluster; }

  protected:
    void do_init(Registry& registry) {
      fjPtrs = registry.accessConst<std::vector<panda::FatJet*>>("fjPtrs");
    }
    void do_execute();
  private:
    std::shared_ptr<const std::vector<panda::FatJet*>> fjPtrs{nullptr};
    std::unique_ptr<fastjet::JetDefinition> jetDef{nullptr};
  };

  // this is not treated as a Op for various reasons
  // the primary being that it does not need access to the tree
  // and that it would make the inheritance complicated because it is
  // intendend to be used by multiple template instances of BaseOp
  class SubRunner {
  public:
    SubRunner(bool ak8, Utils& utils_) : utils(utils_)  {
      double radius = ak8 ? 0.8 : 1.5;
      jetDef.reset(new fastjet::JetDefinition(ak8 ? fastjet::antikt_algorithm :
                                                    fastjet::cambridge_algorithm,
                                              radius));
      tauN.reset(new fastjet::contrib::Njettiness(fastjet::contrib::OnePass_KT_Axes(),
                                                  fastjet::contrib::NormalizedMeasure(1., radius)));
      ecfcalc.reset(new ECFCalculator());

      bool optimalR=true; bool doHTTQ=false;
      double minSJPt=0.; double minCandPt=0.;
      double sjmass=30.; double mucut=0.8;
      double filtR=0.3; int filtN=5;
      int mode=4; double minCandMass=0.;
      double maxCandMass=9999999.; double massRatioWidth=9999999.;
      double minM23Cut=0.; double minM13Cut=0.;
      double maxM13Cut=9999999.;  bool rejectMinR=false;
      htt.reset(new httwrapper::HEPTopTaggerV2(optimalR,doHTTQ,
                                            minSJPt,minCandPt,
                                            sjmass,mucut,
                                            filtR,filtN,
                                            mode,minCandMass,
                                            maxCandMass,massRatioWidth,
                                            minM23Cut,minM13Cut,
                                            maxM13Cut,rejectMinR));
    }

    void run(panda::FatJet& fj);

  private:
    std::unique_ptr<fastjet::JetDefinition> jetDef{nullptr};
    std::unique_ptr<fastjet::contrib::Njettiness>tauN{nullptr};
    std::unique_ptr<ECFCalculator>ecfcalc{nullptr};
    std::unique_ptr<httwrapper::HEPTopTaggerV2>htt{nullptr};
    Utils& utils;
  };

  class FatJetOp : public BaseJetOp {
  public:
    FatJetOp(panda::EventAnalysis& event_,
              Config& cfg_,
              Utils& utils_,
              GeneralTree& gt_,
              int level_=0) :
      BaseJetOp("fatjet", event_, cfg_, utils_, gt_, level_),
      nMaxFJ(analysis.vqqhbb ? 2 : 1),
      fjPtrs(std::v_make_shared<panda::FatJet*>()),
      fatjets(analysis.ak8 ? event.puppiAK8Jets : event.puppiCA15Jets),
      substructure(analysis.recalcECF ? new SubRunner(analysis.ak8, utils) : nullptr) {
        recalcJER = analysis.applyJER || analysis.rerunJER; // if JER is requested, always run it
        recluster = addSubOp<FatJetReclusterOp>();
        jetType = "AK8PFPuppi";
      if      (analysis.year==2016) { // CSVv2 subjet b-tagging in 2016 
        csvL = 0.5426; csvM = 0.8484; 
      } else if (analysis.year==2017) { // DeepCSV subjet b-tagging in 2017
        csvL = 0.1522; csvM = 0.4941;
      }
    }
    virtual ~FatJetOp () { }

    virtual bool on() { return !analysis.genOnly && analysis.fatjet; }

  protected:
    void do_init(Registry& registry) {
      registry.publishConst("fjPtrs", fjPtrs);
      matchLeps = registry.accessConst<std::vector<panda::Lepton*>>("matchLeps");
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
      matchPhos = registry.accessConst<std::vector<panda::Photon*>>("tightPhos");
      dilep = registry.accessConst<TLorentzVector>("dilep");
    }
    void do_reset() {
      fjPtrs->clear();
    }
    void do_execute();
    float getMSDCorr(float,float);
  private:
    void setupJES();

    const int nMaxFJ;
    std::shared_ptr<std::vector<panda::FatJet*>> fjPtrs{nullptr};
    panda::FatJetCollection &fatjets;

    std::shared_ptr<const std::vector<panda::Lepton*>> matchLeps{nullptr};
    std::shared_ptr<const std::vector<panda::Lepton*>> looseLeps{nullptr};
    std::shared_ptr<const std::vector<panda::Photon*>> matchPhos{nullptr};
    std::shared_ptr<const TLorentzVector> dilep{nullptr};
    std::unique_ptr<SubRunner> substructure{nullptr};

    FatJetReclusterOp *recluster{nullptr};
  };

  class FatJetMatchingOp : public AnalysisOp {
  public:
    FatJetMatchingOp(panda::EventAnalysis& event_,
                      Config& cfg_,
                      Utils& utils_,
                      GeneralTree& gt_,
                      int level_=0) :
      AnalysisOp("fjmatching", event_, cfg_, utils_, gt_, level_) { }
    virtual ~FatJetMatchingOp () { }

    virtual bool on() { return !analysis.genOnly && analysis.fatjet && !analysis.isData; }

  protected:
    void do_init(Registry& registry) {
      fjPtrs = registry.accessConst<std::vector<panda::FatJet*>>("fjPtrs");
      if (!analysis.isData)
        genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
    void do_execute();
    void do_reset() { genObjects.clear(); }
  private:
    std::shared_ptr<const std::vector<panda::FatJet*>> fjPtrs{nullptr};
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};

    std::map<const panda::GenParticle*,float> genObjects; // gen particle -> pt
    const panda::GenParticle* matchGen(double eta, double phi, double r2, int pdgid=0) const;
  };

  class HRTagOp : public HROp {
  public:
    HRTagOp(panda::EventAnalysis& event_,
              Config& cfg_,
              Utils& utils_,
              HeavyResTree& gt_,
              int level_=0) :
      HROp("tag", event_, cfg_, utils_, gt_, level_),
      fatjets(analysis.ak8 ? event.puppiAK8Jets : event.puppiCA15Jets),
      substructure(new SubRunner(analysis.ak8, utils)) {
    }
    virtual ~HRTagOp () { }

    virtual bool on() { return true; }

  protected:
    void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
    void do_execute();
  private:

    float getMSDCorr(float,float);
    void fillJet(panda::FatJet&);
    void doSubstructure(panda::FatJet& fj);
    panda::FatJetCollection &fatjets;

    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};
    std::unique_ptr<SubRunner> substructure{nullptr};
  };
}

#endif
