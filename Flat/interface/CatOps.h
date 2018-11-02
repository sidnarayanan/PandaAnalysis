#ifndef CATOPS
#define CATOPS

#include "Operator.h"

namespace pa {
    class CategorizeOp : public AnalysisOp {
      public:
        CategorizeOp(TString name_,
                  panda::EventAnalysis& event_,
                  Config& cfg_,
                  Utils& utils_,
                  GeneralTree& gt_,
                  int level_=0) :
          AnalysisOp(name_+"cat", event_, cfg_, utils_, gt_, level_) { }
        virtual ~CategorizeOp() { }

      protected:
        virtual void do_execute() final { 
          gt.category = categorize();
        }
        virtual int categorize() = 0; 
    };

    class VBFCatOp : public CategorizeOp {
      public:
        VBFCatOp(panda::EventAnalysis& event_,
                  Config& cfg_,
                  Utils& utils_,
                  GeneralTree& gt_,
                  int level_=0) :
          CategorizeOp("vbf", event_, cfg_, utils_, gt_, level_) { }
        virtual ~VBFCatOp() { }

        bool on() override { return analysis.vbf; }
      protected:
        virtual int categorize() override; 
    };
}

#endif
