#ifndef BTAGOPS
#define BTAGOPS

#include "Operator.h"
#include "AnalyzerUtilities.h"

namespace pa {
  class BTagSFOp : public AnalysisOp {
  public:
    BTagSFOp(panda::EventAnalysis& event_, 
              Config& cfg_,                 
              Utils& utils_,                
              GeneralTree& gt_,
              int level_=0) :                 
      AnalysisOp("btagsf", event_, cfg_, utils_, gt_, level_) { }
    virtual ~BTagSFOp () { }

    bool on() { return !analysis.genOnly && analysis.btagSFs && !analysis.isData; }

  protected:
    void do_init(Registry& registry) {
      jesShifts = registry.accessConst<std::vector<JESHandler>>("jesShifts");
    }
    void do_execute();

  private:
    std::shared_ptr<const std::vector<JESHandler>> jesShifts{nullptr};
  };

  class BTagWeightOp : public AnalysisOp {
  public:
    BTagWeightOp(panda::EventAnalysis& event_, 
              Config& cfg_,                 
              Utils& utils_,                
              GeneralTree& gt_,
              int level_=0) :                 
      AnalysisOp("btagweight", event_, cfg_, utils_, gt_, level_) { }
    virtual ~BTagWeightOp () { }

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
