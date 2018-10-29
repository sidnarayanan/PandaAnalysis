#ifndef CATMODS
#define CATMODS

#include "Module.h"

namespace pa {
    class CategorizeMod : public AnalysisMod {
      public:
        CategorizeMod(TString name_,
                  panda::EventAnalysis& event_,
                  Config& cfg_,
                  Utils& utils_,
                  GeneralTree& gt_,
                  int level_=0) :
          AnalysisMod(name_+"cat", event_, cfg_, utils_, gt_, level_) { }
        virtual ~CategorizeMod() { }

      protected:
        virtual void do_execute() final { 
          gt.category = categorize();
        }
        virtual int categorize() = 0; 
    };

    class VBFCatMod : public CategorizeMod {
      public:
        VBFCatMod(panda::EventAnalysis& event_,
                  Config& cfg_,
                  Utils& utils_,
                  GeneralTree& gt_,
                  int level_=0) :
          CategorizeMod("vbf", event_, cfg_, utils_, gt_, level_) { }
        virtual ~VBFCatMod() { }

        bool on() override { return analysis.vbf; }
      protected:
        virtual int categorize() override; 
    };
}

#endif
