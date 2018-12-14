#ifndef JETSOPS
#define JETSOPS

#include "Operator.h"
#include "DeepOps.h"
#include <numeric>

namespace pa {

  class AdJetOp : public AnalysisOp {
  public:
    AdJetOp(panda::EventAnalysis& event_,
	    Config& cfg_,
	    Utils& utils_,
	    GeneralTree& gt_,
	    int level_=0) :
      AnalysisOp("adjet", event_, cfg_, utils_, gt_, level_) { }
    virtual ~AdJetOp() { }
    bool on() { return analysis.hbb; }
  protected:
    void do_init(Registry& registry) {
      currentJES = registry.access<JESHandler*>("currentJES");
    }
    void do_execute(); 
  private:
    std::shared_ptr<JESHandler*> currentJES{nullptr};
  };


  class HbbSystemOp : public AnalysisOp {
  public:
    HbbSystemOp(panda::EventAnalysis& event_,
                 Config& cfg_,
                 Utils& utils_,
                 GeneralTree& gt_,
                 int level_=0) :
      AnalysisOp("hbbsystem", event_, cfg_, utils_, gt_, level_),
      hbbdJet(std::make_shared<JetWrapper*>(nullptr)) {
        deepreg = addSubOp<BRegDeepOp>();
        bdtreg = addSubOp<BRegBDTOp>();
      }
    virtual ~HbbSystemOp() { }

    bool on() { return analysis.hbb; }
  protected:
    void do_init(Registry& registry) {
      currentJES = registry.access<JESHandler*>("currentJES");
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
      dilep = registry.accessConst<TLorentzVector>("dilep");
      registry.publish("higgsDaughterJet", hbbdJet);
    }
    void do_execute();
    void do_reset() { btagsorted.clear(); }
  private:
    std::shared_ptr<JESHandler*> currentJES{nullptr};
    std::shared_ptr<const std::vector<panda::Lepton*>> looseLeps{nullptr};
    std::shared_ptr<const TLorentzVector> dilep{nullptr};

    std::vector<JetWrapper*> btagsorted;

    std::shared_ptr<JetWrapper*> hbbdJet;
    BRegDeepOp *deepreg{nullptr};
    BRegBDTOp *bdtreg{nullptr};
  };

  class JetFlavorOp : public AnalysisOp {
  public:
    JetFlavorOp(panda::EventAnalysis& event_,
                 Config& cfg_,
                 Utils& utils_,
                 GeneralTree& gt_,
                 int level_=0) :
      AnalysisOp("jetflavor", event_, cfg_, utils_, gt_, level_) { }
    virtual ~JetFlavorOp () {}

    bool on() { return !analysis.isData && (analysis.jetFlavorPartons || analysis.jetFlavorJets); }
  protected:
    void do_init(Registry& registry) {
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
      if (!analysis.isData)
        genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
    void do_execute();
  private:
    std::shared_ptr<std::vector<JESHandler>> jesShifts{nullptr};
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};

    void partonFlavor(JetWrapper&);
    void clusteredFlavor(JetWrapper&);
  };

  class IsoJetOp : public AnalysisOp {
  public:
    IsoJetOp(panda::EventAnalysis& event_,
              Config& cfg_,
              Utils& utils_,
              GeneralTree& gt_,
              int level_=0) :
      AnalysisOp("isojet", event_, cfg_, utils_, gt_, level_) { }
    virtual ~IsoJetOp () {}

    bool on() { return analysis.fatjet; }
  protected:
    void do_init(Registry& registry) {
      currentJet = registry.access<JetWrapper*>("currentJet");
      currentJES = registry.access<JESHandler*>("currentJES");

      fjPtrs = registry.accessConst<std::vector<panda::FatJet*>>("fjPtrs");
    }
    void do_execute();
  private:
    std::shared_ptr<JetWrapper*> currentJet{nullptr}; // shared ptr to a bare address
    std::shared_ptr<JESHandler*> currentJES{nullptr};
    std::shared_ptr<const std::vector<panda::FatJet*>> fjPtrs{nullptr};
  };

  class BJetRegOp : public AnalysisOp {
  public:
    BJetRegOp(panda::EventAnalysis& event_,
                  Config& cfg_,
                  Utils& utils_,
                  GeneralTree& gt_,
                  int level_=0) :
      AnalysisOp("bjetreg", event_, cfg_, utils_, gt_, level_) { }
    virtual ~BJetRegOp () {}

    bool on() { return analysis.bjetBDTReg || analysis.bjetDeepReg || analysis.bjetRegTraining; }
  protected:
    void do_init(Registry& registry) {
      currentJet = registry.access<JetWrapper*>("currentJet");
      currentJES = registry.access<JESHandler*>("currentJES");
    }
    void do_execute();
    void do_reset() { energies.clear(); }
  private:
    struct Energies {
      static const std::vector<double> dr2_bins;
      enum pftype {
        pem, pch, pmu, pne, pN
      };
      using contents = std::array<std::vector<TLorentzVector>, static_cast<int>(shiftjetrings::N)>;
      std::array<contents, pN> pf; // pf[pftype][bin_idx]
      float jet_e; 

      void clear() {
        for (auto& i : pf)
          for (auto& j : i)
            j.clear();
      }
      float get_e(int bin, pftype pft) {
        float sum = 0;
        for (auto& v : pf[pft][bin])
          sum += v.E();
        return sum / jet_e;
      }

      typedef double (TLorentzVector::*vec_func) () const;
      template <long unsigned int I>
      void get_moments(int pft, vec_func f, std::array<float,I>& moments, float shift=0) {
        std::fill(moments.begin(), moments.end(), 0);
        // first get the mean
        float sumw{0};
        std::vector<std::pair<float,float>> x;
        for (auto& bin : pf[pft]) {
          x.reserve(x.size() + bin.size());
          for (auto& v : bin) {
            x.emplace_back((v.*f)() - shift, v.E());
            sumw += v.E();
          }
        }

        if (sumw == 0)
          return;

        for (auto& xx : x)
          moments[0] += xx.first * xx.second; 
        moments[0] /= sumw;

        for (unsigned ex = 1; ex != I; ++ ex) {
          for (auto& xx : x) {
            // dividing by sumw here is not optimal but prevents a very large sum
            moments[ex] += std::pow(xx.second / sumw * (xx.first - moments[0]), ex + 1);
          }
          moments[ex] = std::pow(std::fabs(moments[ex]), 1./(ex + 1));
        }
      }
    };

    Energies energies;
    std::shared_ptr<JetWrapper*> currentJet{nullptr}; // shared ptr to a bare address
    std::shared_ptr<JESHandler*> currentJES{nullptr};
  };

  class VBFSystemOp : public AnalysisOp {
  public:
    VBFSystemOp(panda::EventAnalysis& event_,
                 Config& cfg_,
                 Utils& utils_,
                 GeneralTree& gt_,
                 int level_=0) :
      AnalysisOp("vbfsystem", event_, cfg_, utils_, gt_, level_) { }
    virtual ~VBFSystemOp () {}

    bool on() { return analysis.vbf; }
  protected:
    void do_init(Registry& registry) {
      currentJES = registry.access<JESHandler*>("currentJES");
      fjPtrs = registry.accessConst<std::vector<panda::FatJet*>>("fjPtrs", true); 
    }
    void do_execute();
  private:
    std::shared_ptr<JESHandler*> currentJES{nullptr};
    std::shared_ptr<const std::vector<panda::FatJet *>> fjPtrs{nullptr};
  };

  class BaseJetOp : public AnalysisOp {
  public:
    BaseJetOp(TString name,
               panda::EventAnalysis& event_,
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) :
      AnalysisOp(name, event_, cfg_, utils_, gt_, level_),
      recalcJER(false) {
        if (analysis.year == 2016) {
          jecV = "V4"; jecReco = "23Sep2016";
          campaign = "Summer16";
          jerV = "Spring16_25nsV10";
          eraGroups = {"BCD","EF","G","H"};
          spacer = "";
          if (analysis.useDeepCSV) { 
            csvL = 0.2219; csvM = 0.6324; 
          } else { 
            csvL = 0.5426; csvM = 0.8484; 
          }
        } 
	else if (analysis.year == 2017 || analysis.year == 2018) {
          jecV = "V8"; jecReco = "17Nov2017";
          //jecV = "V32"; jecReco = "17Nov2017"; // missing V32 for now
          campaign = "Fall17";
          jerV = "Fall17_25nsV1";
          eraGroups = {"B","C","D","E","F"};
          // eraGroups = {"B","C","DE","F"}; // missing V32 for now
          spacer = "_";
          if (analysis.useDeepCSV) { 
            csvL = 0.1522; csvM = 0.4941; 
          } else { 
            csvL = 0.5803; csvM = 0.8838; 
          }
        }
      }
    virtual ~BaseJetOp() { }
    bool csvLoose(float csv) { return csv > csvL; }
    bool csvMed(float csv) { return csv > csvM; }
  protected:
    bool recalcJER;
    virtual void do_execute() = 0;
    virtual void do_readData(TString path);
    JetWrapper shiftJet(const panda::Jet& jet, shiftjes shift, bool smear=false);

    std::map<TString,std::unique_ptr<FactorizedJetCorrector>> scales; // era/MC -> scale
    std::map<TString,std::vector<std::shared_ptr<JetCorrectionUncertainty>>> scaleUncs; // era/MC -> (src -> unc)
    std::unique_ptr<JERReader> jer{nullptr}; //!< fatjet jet energy resolution reader

    std::vector<std::shared_ptr<JetCorrectionUncertainty>> *scaleUnc  {nullptr}; // src -> unc
    FactorizedJetCorrector *scale{nullptr};

    TString jecV, jecReco, jetType, campaign, spacer, jerV;
    std::vector<TString> eraGroups;
    float csvL, csvM;
  private:
    void setScaleUnc(TString, TString);
  };

  class JetOp : public BaseJetOp {
  public:
    JetOp(panda::EventAnalysis& event_,
           Config& cfg_,
           Utils& utils_,
           GeneralTree& gt_,
           int level_=0) :
      BaseJetOp("jet", event_, cfg_, utils_, gt_, level_),
      currentJet(std::make_shared<JetWrapper*>(nullptr)),
      currentJES(std::make_shared<JESHandler*>(nullptr)) {
        ak4Jets = &(event.chsAK4Jets);
        recalcJER = analysis.rerunJER; 

        isojet = addSubOp<IsoJetOp>();
        bjetreg = addSubOp<BJetRegOp>();
        vbf = addSubOp<VBFSystemOp>();
        hbb = addSubOp<HbbSystemOp>();
	adjet = addSubOp<AdJetOp>();

        jetType = "AK4PFchs";
    }
    virtual ~JetOp () { }

    virtual bool on() { return !analysis.genOnly; }

  protected:
    void do_init(Registry& registry) {
      registry.publish("currentJet", currentJet);
      registry.publish("currentJES", currentJES);
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
      matchLeps = registry.accessConst<std::vector<panda::Lepton*>>("matchLeps");
      matchPhos = registry.accessConst<std::vector<panda::Photon*>>("tightPhos");
      matchVeryLoosePhos = registry.accessConst<std::vector<panda::Photon*>>("veryLoosePhos");
    }
    void do_execute();

  private:
    IsoJetOp *isojet{nullptr};
    BJetRegOp *bjetreg{nullptr};
    VBFSystemOp *vbf{nullptr};
    HbbSystemOp *hbb{nullptr};
    AdJetOp *adjet{nullptr};

    std::shared_ptr<std::vector<JESHandler>> jesShifts{nullptr};

    std::shared_ptr<const std::vector<panda::Lepton*>> matchLeps{nullptr};
    std::shared_ptr<const std::vector<panda::Photon*>> matchPhos{nullptr};
    std::shared_ptr<const std::vector<panda::Photon*>> matchVeryLoosePhos{nullptr};

    panda::JetCollection *ak4Jets{nullptr};

    std::shared_ptr<JetWrapper*> currentJet; // shared ptr to a bare address
    std::shared_ptr<JESHandler*> currentJES;

    void setupJES();
    void varyJES();
  };
}

#endif
