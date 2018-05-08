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
    ~HbbSystemMod () {
      delete bjetregReader;
      delete[] bjetreg_vars;
    }

    bool on() { return analysis.hbb; }
  protected:
    void do_readData(TString dirPath);
    void do_init(Registry& registry) {
      currentJES = registry.access<JESHandler*>("currentJES");
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps");
      registry.publishConst("btagsortedjets", &btagsorted);
    }
    void do_execute();
  private:
    JESHandler **currentJES{nullptr};
    std::vector<const JetWrapper*> btagsorted;
    const std::vector<panda::Lepton*>* looseLeps{nullptr};

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
    ~JetFlavorMod () {}

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
    ~IsoJetMod () {}

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
    ~BJetRegMod () {}

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
    ~VBFSystemMod () {}

    bool on() { return analysis.vbf; }
  protected:
    void do_init(Registry& registry) {
      currentJES = registry.access<JESHandler*>("currentJES");
    }
    void do_execute();
  private:
    JESHandler **currentJES{nullptr};
  };


  class JetMod : public AnalysisMod {
  public: 
    JetMod(panda::EventAnalysis& event_, 
           Config& cfg_,                 
           Utils& utils_,                
           GeneralTree& gt_) :                 
      AnalysisMod("jet", event_, cfg_, utils_, gt_) { 
        ak4Jets = &(event.chsAK4Jets); 

        flavor = new JetFlavorMod(event_, cfg_, utils_, gt_); subMods.push_back(flavor);
        isojet = new IsoJetMod(event_, cfg_, utils_, gt_); subMods.push_back(isojet);
        bjetreg = new BJetRegMod(event_, cfg_, utils_, gt_); subMods.push_back(bjetreg);
        vbf = new VBFSystemMod(event_, cfg_, utils_, gt_); subMods.push_back(vbf);
        hbb = new HbbSystemMod(event_, cfg_, utils_, gt_); subMods.push_back(hbb);
      }
    ~JetMod () { 
      delete ak4JERReader;
      for (auto& iter : ak4ScaleReader)
        delete iter.second;
      for (auto& iter : ak4UncReader)
        delete iter.second;
    }

    virtual bool on() { return !analysis.genOnly; }
    
  protected:
    void do_readData(TString path);
    void do_init(Registry& registry) {
      registry.publish("currentJet", &currentJet);
      registry.publish("currentJES", &currentJES);
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
      matchLeps = registry.accessConst<std::vector<panda::Lepton*>>("matchLeps");
      looseLeps = registry.accessConst<std::vector<panda::Lepton*>>("looseLeps"); // needed for ZH spin
      matchPhos = registry.accessConst<std::vector<panda::Photon*>>("tightPhos");
      dilep = registry.accessConst<TLorentzVector*>("dilep"); // needed for ZH spin
    }
    void do_execute();  

  private:
    JetFlavorMod *flavor{nullptr};
    IsoJetMod *isojet{nullptr};
    BJetRegMod *bjetreg{nullptr};
    VBFSystemMod *vbf{nullptr};
    HbbSystemMod *hbb{nullptr};

    std::map<TString,FactorizedJetCorrector*> ak4ScaleReader; //!< calculate JES on the fly
    std::map<TString,JetCorrectionUncertainty*> ak4UncReader; //!< calculate JES unc on the fly
    JERReader *ak4JERReader{nullptr}; //!< fatjet jet energy resolution reader
    JetCorrectionUncertainty *uncReaderAK4  {nullptr};        
    FactorizedJetCorrector   *scaleReaderAK4{nullptr};        
    std::vector<JESHandler>* jesShifts{nullptr}; 

    const std::vector<panda::Lepton*>* matchLeps{nullptr};
    const std::vector<panda::Lepton*>* looseLeps{nullptr};
    const std::vector<panda::Photon*>* matchPhos{nullptr};
    const TLorentzVector *dilep{nullptr};

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
    ~SoftActivityMod () { 
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
    ~GenJetNuMod () { delete jetDef; }

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
