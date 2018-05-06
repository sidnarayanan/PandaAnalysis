#ifndef BTAGMODS
#define BTAGMODS

#include "Module.h"
#include "AnalyzerUtilities.h"

namespace pa {
  class BTagSFMod : public AnalysisMod {
  public:
    BTagSFMod(panda::EventAnalysis& event_, 
              Config& cfg_,                 
              Utils& utils_,                
              GeneralTree& gt_) :                 
      AnalysisMod("btagsf", event_, cfg_, utils_, gt_) { }
    ~BTagSFMod () { }

    bool on() { return !analysis.genOnly && analysis.btagSFs; }

  protected:
    void do_init(Registry& registry) {
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
    }
    void do_execute();

  private:
    const std::vector<JESHandler> *jesShifts{nullptr};
  };

  class BTagWeightMod : public AnalysisMod {
  public:
    BTagWeightMod(panda::EventAnalysis& event_, 
              Config& cfg_,                 
              Utils& utils_,                
              GeneralTree& gt_) :                 
      AnalysisMod("btagweight", event_, cfg_, utils_, gt_) { }
    ~BTagWeightMod () { }

    bool on() { return !analysis.genOnly && analysis.btagWeights; }

  protected:
    void do_init(Registry& registry) {
      jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
    }
    void do_execute();

  private:
    const std::vector<JESHandler> *jesShifts{nullptr};
  };
 
}

#endif
