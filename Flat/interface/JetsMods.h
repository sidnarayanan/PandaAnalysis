#ifndef JETSMODS
#define JETSMODS

#include "Module.h"
#include "DeepMods.h"
#include <numeric>

namespace pa {
  class HbbSystemMod : public AnalysisMod {
  public:
    HbbSystemMod(panda::EventAnalysis& event_,
                 Config& cfg_,
                 Utils& utils_,
                 GeneralTree& gt_,
                 int level_=0) :
      AnalysisMod("hbbsystem", event_, cfg_, utils_, gt_, level_),
      hbbdJet(std::make_shared<JetWrapper*>(nullptr)) {
        deepreg = addSubMod<BRegDeepMod>();
        bdtreg = addSubMod<BRegBDTMod>();
      }
    virtual ~HbbSystemMod() { }

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
    BRegDeepMod *deepreg{nullptr};
    BRegBDTMod *bdtreg{nullptr};
  };

  class JetFlavorMod : public AnalysisMod {
  public:
    JetFlavorMod(panda::EventAnalysis& event_,
                 Config& cfg_,
                 Utils& utils_,
                 GeneralTree& gt_,
                 int level_=0) :
      AnalysisMod("jetflavor", event_, cfg_, utils_, gt_, level_) { }
    virtual ~JetFlavorMod () {}

    bool on() { return analysis.jetFlavorPartons || analysis.jetFlavorJets; }
  protected:
    void do_init(Registry& registry) {
      currentJet = registry.access<JetWrapper*>("currentJet");
      if (!analysis.isData)
        genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
    void do_execute();
  private:
    std::shared_ptr<JetWrapper*> currentJet{nullptr}; // shared ptr to a bare address
    std::shared_ptr<const std::vector<panda::Particle*>> genP{nullptr};

    void partonFlavor();
    void clusteredFlavor();
  };

  class IsoJetMod : public AnalysisMod {
  public:
    IsoJetMod(panda::EventAnalysis& event_,
              Config& cfg_,
              Utils& utils_,
              GeneralTree& gt_,
              int level_=0) :
      AnalysisMod("isojet", event_, cfg_, utils_, gt_, level_) { }
    virtual ~IsoJetMod () {}

    bool on() { return analysis.fatjet; }
  protected:
    void do_init(Registry& registry) {
      currentJet = registry.access<JetWrapper*>("currentJet");
      currentJES = registry.access<JESHandler*>("currentJES");

      fj1 = registry.accessConst<const panda::FatJet*>("fj1");
    }
    void do_execute();
  private:
    std::shared_ptr<JetWrapper*> currentJet{nullptr}; // shared ptr to a bare address
    std::shared_ptr<JESHandler*> currentJES{nullptr};
    std::shared_ptr<const panda::FatJet* const> fj1{nullptr};
  };

  class BJetRegMod : public AnalysisMod {
  public:
    BJetRegMod(panda::EventAnalysis& event_,
                  Config& cfg_,
                  Utils& utils_,
                  GeneralTree& gt_,
                  int level_=0) :
      AnalysisMod("bjetreg", event_, cfg_, utils_, gt_, level_) { }
    virtual ~BJetRegMod () {}

    bool on() { return analysis.bjetBDTReg || analysis.bjetDeepReg; }
  protected:
    void do_init(Registry& registry) {
      currentJet = registry.access<JetWrapper*>("currentJet");
      currentJES = registry.access<JESHandler*>("currentJES");
    }
    void do_execute();
    void do_reset() { energies.clear(); }
  private:
    struct Energies {
      static const std::vector<double> dr_bins;
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
      void get_moments(int pft, vec_func f, std::array<float,I>& moments ) {
        std::fill(moments.begin(), moments.end(), 0);
        // first get the mean
        float sumw{0};
        std::vector<std::pair<float,float>> x;
        for (auto& bin : pf[pft]) {
          x.reserve(x.size() + bin.size());
          for (auto& v : bin) {
            x.emplace_back((v.*f)(),v.E());
            sumw += v.E();
          }
        }

        if (sumw == 0)
          return;

        for (auto& xx : x)
          moments[0] += (xx.first * xx.second); 
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

  class VBFSystemMod : public AnalysisMod {
  public:
    VBFSystemMod(panda::EventAnalysis& event_,
                 Config& cfg_,
                 Utils& utils_,
                 GeneralTree& gt_,
                 int level_=0) :
      AnalysisMod("vbfsystem", event_, cfg_, utils_, gt_, level_) { }
    virtual ~VBFSystemMod () {}

    bool on() { return analysis.vbf; }
  protected:
    void do_init(Registry& registry) {
      currentJES = registry.access<JESHandler*>("currentJES");
      fj1 = registry.accessConst<const panda::FatJet*>("fj1", true); 
    }
    void do_execute();
  private:
    std::shared_ptr<JESHandler*> currentJES{nullptr};
    std::shared_ptr<const panda::FatJet *const> fj1{nullptr};
  };

  class BaseJetMod : public AnalysisMod {
  public:
    BaseJetMod(TString name,
               panda::EventAnalysis& event_,
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_,
               int level_=0) :
      AnalysisMod(name, event_, cfg_, utils_, gt_, level_) {
        if (analysis.year == 2016) {
          jecV = "V4"; jecReco = "23Sep2016";
          campaign = "Summer16";
          jerV = "Spring16_25nsV10";
          eraGroups = {"BCD","EF","G","H"};
          spacer = "";
          csvL = 0.5426; csvM = 0.8484;
        } else {
          jecV = "V8"; jecReco = "17Nov2017";
          campaign = "Fall17";
          jerV = "Fall17_25nsV1";
          eraGroups = {"B","C","D","E","F"};
          spacer = "_";
          csvL = 0.2219; csvM = 0.6324;
        }
      }
    virtual ~BaseJetMod () { }
    bool csvLoose (float csv) { return csv > csvL; }
    bool csvMed (float csv) { return csv > csvM; }
  protected:
    virtual void do_execute() = 0;
    virtual void do_readData(TString path);
    JetWrapper shiftJet(const panda::Jet& jet, shiftjes shift, bool smear=false, bool recalcSmear=false);

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

  class JetMod : public BaseJetMod {
  public:
    JetMod(panda::EventAnalysis& event_,
           Config& cfg_,
           Utils& utils_,
           GeneralTree& gt_,
           int level_=0) :
      BaseJetMod("jet", event_, cfg_, utils_, gt_, level_),
      currentJet(std::make_shared<JetWrapper*>(nullptr)),
      currentJES(std::make_shared<JESHandler*>(nullptr)) {
        ak4Jets = &(event.chsAK4Jets);

        flavor = addSubMod<JetFlavorMod>();
        isojet = addSubMod<IsoJetMod>();
        bjetreg = addSubMod<BJetRegMod>();
        vbf = addSubMod<VBFSystemMod>();
        hbb = addSubMod<HbbSystemMod>();

        jetType = "AK4PFchs";
    }
    virtual ~JetMod () { }

    virtual bool on() { return !analysis.genOnly; }

  protected:
    void do_init(Registry& registry) {
      registry.publish("currentJet", currentJet);
      registry.publish("currentJES", currentJES);
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
      matchLeps = registry.accessConst<std::vector<panda::Lepton*>>("matchLeps");
      matchPhos = registry.accessConst<std::vector<panda::Photon*>>("tightPhos");
    }
    void do_execute();

  private:
    JetFlavorMod *flavor{nullptr};
    IsoJetMod *isojet{nullptr};
    BJetRegMod *bjetreg{nullptr};
    VBFSystemMod *vbf{nullptr};
    HbbSystemMod *hbb{nullptr};

    std::shared_ptr<std::vector<JESHandler>> jesShifts{nullptr};

    std::shared_ptr<const std::vector<panda::Lepton*>> matchLeps{nullptr};
    std::shared_ptr<const std::vector<panda::Photon*>> matchPhos{nullptr};

    panda::JetCollection *ak4Jets{nullptr};

    std::shared_ptr<JetWrapper*> currentJet; // shared ptr to a bare address
    std::shared_ptr<JESHandler*> currentJES;

    void setupJES();
    void varyJES();
  };
}

#endif
