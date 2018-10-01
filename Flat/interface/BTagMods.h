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
              GeneralTree& gt_,
              int level_=0) :                 
      AnalysisMod("btagsf", event_, cfg_, utils_, gt_, level_) { }
    virtual ~BTagSFMod () { }

    bool on() { return !analysis.genOnly && analysis.btagSFs && !analysis.isData; }

  protected:
    void do_init(Registry& registry) {
      jesShifts = registry.accessConst<std::vector<JESHandler>>("jesShifts");
    }
    void do_execute();

  private:
    std::shared_ptr<const std::vector<JESHandler>> jesShifts{nullptr};
  };

  class BTagWeightMod : public AnalysisMod {
  public:
    BTagWeightMod(panda::EventAnalysis& event_, 
              Config& cfg_,                 
              Utils& utils_,                
              GeneralTree& gt_,
              int level_=0) :                 
      AnalysisMod("btagweight", event_, cfg_, utils_, gt_, level_) { }
    virtual ~BTagWeightMod () { }

    bool on() { return !analysis.genOnly && analysis.btagWeights && !analysis.isData; }

  protected:
    void do_init(Registry& registry) {
      jesShifts = registry.accessConst<std::vector<JESHandler>>("jesShifts");
    }
    void do_execute();

  private:
    std::shared_ptr<const std::vector<JESHandler>> jesShifts{nullptr};
  };
 
}

#endif
