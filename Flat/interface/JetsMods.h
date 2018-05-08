#ifndef JETSMODS
#define JETSMODS

#include "Module.h"
#include "TMVA/Reader.h"

namespace pa {

  class HbbSystemMod : public AnalysisMod {
  public:
    HbbSystemMod(panda::EventAnalysis& event_, 
                 Config& cfg_,                 
                 Utils& utils_,                
                 GeneralTree& gt_) :                 
      AnalysisMod("hbbsystem", event_, cfg_, utils_, gt_) { }
    virtual ~HbbSystemMod () {
      delete bjetregReader;
      delete[] bjetreg_vars;
    }

    bool on() { return analysis.hbb; }
  protected:
    void do_readData(TString dirPath);
    void do_init(Registry& registry) {
      currentJES = registry.access<JESHandler*>("currentJES");
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps"); // needed for W reco and ZH spin
      registry.publishConst("btagsortedjets", &btagsorted);
      dilep = registry.accessConst<TLorentzVector>("dilep"); // needed for ZH spin
    }
    void do_execute();
  private:
    JESHandler **currentJES{nullptr};
    std::vector<const JetWrapper*> btagsorted;
    const std::vector<panda::Lepton*>* looseLeps{nullptr};
    const TLorentzVector *dilep{nullptr};

    TMVA::Reader *bjetregReader{nullptr}; 
    float *bjetreg_vars{nullptr};
  };

  class JetFlavorMod : public AnalysisMod {
  public:
    JetFlavorMod(panda::EventAnalysis& event_, 
                 Config& cfg_,                 
                 Utils& utils_,                
                 GeneralTree& gt_) :                 
      AnalysisMod("jetflavor", event_, cfg_, utils_, gt_) { }
    virtual ~JetFlavorMod () {}

    bool on() { return analysis.jetFlavorPartons || analysis.jetFlavorJets; }
  protected:
    void do_init(Registry& registry) {
      currentJet = registry.access<JetWrapper*>("currentJet");
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
    void do_execute();
  private:
    JetWrapper **currentJet{nullptr};
    const std::vector<panda::Particle*> *genP;

    void partonFlavor();
    void clusteredFlavor();
  };

  class IsoJetMod : public AnalysisMod {
  public:
    IsoJetMod(panda::EventAnalysis& event_, 
              Config& cfg_,                 
              Utils& utils_,                
              GeneralTree& gt_) :                 
      AnalysisMod("isojet", event_, cfg_, utils_, gt_) { }
    virtual ~IsoJetMod () {}

    bool on() { return analysis.fatjet; }
  protected:
    void do_init(Registry& registry) {
      currentJet = registry.access<JetWrapper*>("currentJet");
      currentJES = registry.access<JESHandler*>("currentJES");
      fj1 = registry.accessConst<panda::FatJet*>("fj1");
    }
    void do_execute();
  private:
    JetWrapper **currentJet{nullptr};
    JESHandler **currentJES{nullptr};
    panda::FatJet *const *fj1{nullptr}; 
  };

  class BJetRegMod : public AnalysisMod {
  public:
    BJetRegMod(panda::EventAnalysis& event_, 
                  Config& cfg_,                 
                  Utils& utils_,                
                  GeneralTree& gt_) :                 
      AnalysisMod("bjetreg", event_, cfg_, utils_, gt_) { }
    virtual ~BJetRegMod () {}

    bool on() { return analysis.bjetRegression; }
  protected:
    void do_init(Registry& registry) {
      currentJet = registry.access<JetWrapper*>("currentJet");
      currentJES = registry.access<JESHandler*>("currentJES");
    }
    void do_execute();
  private:
    JetWrapper **currentJet{nullptr};
    JESHandler **currentJES{nullptr};
  };

  class VBFSystemMod : public AnalysisMod {
  public:
    VBFSystemMod(panda::EventAnalysis& event_, 
                 Config& cfg_,                 
                 Utils& utils_,                
                 GeneralTree& gt_) :                 
      AnalysisMod("vbfsystem", event_, cfg_, utils_, gt_) { }
    virtual ~VBFSystemMod () {}

    bool on() { return analysis.vbf; }
  protected:
    void do_init(Registry& registry) {
      currentJES = registry.access<JESHandler*>("currentJES");
    }
    void do_execute();
  private:
    JESHandler **currentJES{nullptr};
  };

  class BaseJetMod : public AnalysisMod {
  public: 
    BaseJetMod(TString name,
               panda::EventAnalysis& event_,
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_) :
      AnalysisMod(name, event_, cfg_, utils_, gt_) { 
        if (analysis.year == 2016) {
          jecV = "V4"; jecReco = "23Sep2016"; 
          campaign = "Summer16";
          jerV = "Spring16_25nsV10";
          eraGroups = {"BCD","EF","G","H"};
          spacer = "";
        } else {
          jecV = "V8"; jecReco = "17Nov2017"; 
          campaign = "Fall17";
          jerV = "Fall17_25nsV1";
          eraGroups = {"B","C","D","E","F"};
          spacer = "_";
        }
      }
    virtual ~BaseJetMod () { 
      delete jer;
      for (auto& iter : scales)
        delete iter.second;
      for (auto& iter : scaleUncs)
        for (auto* p : iter.second)
          delete p;
    }
  protected:
    virtual void do_execute() = 0;
    virtual void do_readData(TString path);
    JetWrapper shiftJet(const panda::Jet& jet, shiftjes shift, bool smear=false);

    std::map<TString,FactorizedJetCorrector*> scales; // era/MC -> scale 
    std::map<TString,std::vector<JetCorrectionUncertainty*>> scaleUncs; // era/MC -> (src -> unc)
    JERReader *jer{nullptr}; //!< fatjet jet energy resolution reader
    std::vector<JetCorrectionUncertainty*> *scaleUnc  {nullptr}; // src -> unc 
    FactorizedJetCorrector   *scale{nullptr};        
    
    TString jecV, jecReco, jetType, campaign, spacer, jerV;
    std::vector<TString> eraGroups;
  private:
    void setScaleUnc(TString, TString);
  };

  class JetMod : public BaseJetMod {
  public: 
    JetMod(panda::EventAnalysis& event_, 
           Config& cfg_,                 
           Utils& utils_,                
           GeneralTree& gt_) :                 
      BaseJetMod("jet", event_, cfg_, utils_, gt_) { 
        ak4Jets = &(event.chsAK4Jets); 

        flavor = new JetFlavorMod(event_, cfg_, utils_, gt_); subMods.push_back(flavor);
        isojet = new IsoJetMod(event_, cfg_, utils_, gt_); subMods.push_back(isojet);
        bjetreg = new BJetRegMod(event_, cfg_, utils_, gt_); subMods.push_back(bjetreg);
        vbf = new VBFSystemMod(event_, cfg_, utils_, gt_); subMods.push_back(vbf);
        hbb = new HbbSystemMod(event_, cfg_, utils_, gt_); subMods.push_back(hbb);

        jetType = "AK4PFchs";
      }
    virtual ~JetMod () { }

    virtual bool on() { return !analysis.genOnly; }
    
  protected:
    void do_init(Registry& registry) {
      registry.publish("currentJet", &currentJet);
      registry.publish("currentJES", &currentJES);
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

    std::vector<JESHandler>* jesShifts{nullptr}; 

    const std::vector<panda::Lepton*>* matchLeps{nullptr};
    const std::vector<panda::Photon*>* matchPhos{nullptr};

    panda::JetCollection *ak4Jets{nullptr};

    JetWrapper *currentJet{nullptr};
    JESHandler *currentJES{nullptr};

    void setupJES();
    void varyJES();
  };


  class SoftActivityMod : public AnalysisMod {
  public: 
    SoftActivityMod(panda::EventAnalysis& event_, 
                    Config& cfg_,                 
                    Utils& utils_,                
                    GeneralTree& gt_) :                 
      AnalysisMod("softactivity", event_, cfg_, utils_, gt_) { 
      }
    virtual ~SoftActivityMod () { 
      delete jetDefSoftTrack;
    }

    virtual bool on() { return !analysis.genOnly && analysis.hbb; }
    
  protected:
    void do_init(Registry& registry) {
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
      jetDefSoftTrack = new fastjet::JetDefinition(fastjet::antikt_algorithm,0.4);
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
    }
    void do_execute();  

  private:
    std::vector<JESHandler>* jesShifts{nullptr}; 
    fastjet::JetDefinition* jetDefSoftTrack       {nullptr};
    const std::vector<panda::Lepton*>* looseLeps{nullptr};
  };


  class GenJetNuMod : public AnalysisMod {
  public:
    GenJetNuMod(panda::EventAnalysis& event_, 
                 Config& cfg_,                 
                 Utils& utils_,                
                 GeneralTree& gt_) :                 
      AnalysisMod("genjetnu", event_, cfg_, utils_, gt_) { 
        jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm,0.4);
      }
    virtual ~GenJetNuMod () { delete jetDef; }

    bool on() { return !analysis.isData && analysis.reclusterGen && analysis.hbb; }
  protected:
    void do_init(Registry& registry) {
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
    }
    void do_execute();
  private:
    std::vector<JESHandler>* jesShifts{nullptr}; 
    const std::vector<panda::Particle*> *genP;
    fastjet::JetDefinition *jetDef{nullptr}; 
  };
}

#endif
