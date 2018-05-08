#ifndef HBBMODS
#define HBBMODS

#include "Module.h"

namespace pa {
  class HbbMiscMod : public AnalysisMod {
  public: 
    HbbMiscMod(panda::EventAnalysis& event_, 
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_) : 
      AnalysisMod("hbbmisc", event_, cfg_, utils_, gt_) { }
    virtual ~HbbMiscMod () { }

    virtual bool on() { return !analysis.genOnly && analysis.hbb; }
    
  protected:
    void do_execute(); 
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
