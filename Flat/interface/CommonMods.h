#include "Module.cc"
#include "AnalyzerUtilities.h"

namespace pa {
  class TriggerMod : public AnalysisMod {
  public:
    TriggerMod(const panda::EventAnalysis& event_, 
               const Config& cfg_,                 
               const Utils& utils_,                
               GeneralTree& gt_) :                 
      AnalysisMod("trigger", event_, cfg_, utils_, gt_) { }
    ~TriggerMod () { }

  protected:
    void do_initalize(Registry& registry);
    void do_execute();

  private:
    auto triggerHandlers = std::vector<TriggerHandler>(kNTrig) 
  };

  class GlobalMod : public AnalysisMod {
  public:
    GlobalMod(const panda::EventAnalysis& event_, 
               const Config& cfg_,                 
               const Utils& utils_,                
               GeneralTree& gt_) :                 
      AnalysisMod("global", event_, cfg_, utils_, gt_) { }
    ~GlobalMod () { }
    
  protected:
    void do_init(Registry& registry);
    void do_execute();  
  private:
    auto jesShifts = std::vector<JESHandler>(jes2i(shiftjes::N)); 
  };

  class TemplateDupModMod : public AnalysisMod {
  public:
    TemplateDupModMod(const panda::EventAnalysis& event_, 
                      const Config& cfg_,                 
                      const Utils& utils_,                
                      GeneralTree& gt_) :                 
      AnalysisMod("gendup", event_, cfg_, utils_, gt_) { }
    ~TemplateDupModMod () { 
      if (owned) {
        delete fsGenP;
        delete nfsGenP;
      }
    }
    
  protected:
    void do_init(Registry& registry);
    void do_execute();  
  private:
    std::vector<panda::Particle*> *fsGenP{nullptr}, *nfsGenP{nullptr};
    bool owned{true};
  };
}
