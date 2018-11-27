#ifndef HRAnalyzer_h
#define HRAnalyzer_h

#include "AnalyzerUtilities.h"
#include "Analyzer.h"
#include "HeavyResTree.h"
#include "Selection.h"
#include "Operator.h"
#include "PandaAnalysis/Flat/interface/Common.h"

#include "CommonOps.h"
#include "FatJetsOps.h"


namespace pa {

    class HRAnalyzer : public Analyzer<HeavyResTree> {
    public :
        HRAnalyzer(Analysis* a, int debug_=0);
        ~HRAnalyzer();
        void Reset();
        void Run();
        void Terminate();

    private:

        //////////////////////////////////////////////////////////////////////////////////////
        std::unique_ptr<HRGenPOp> gen{nullptr};
        std::unique_ptr<FatJetReclusterOp<HeavyResTree>> fj{nullptr};
        std::unique_ptr<HRTagOp> op{nullptr};
        Config cfg;
        Utils utils; 

        panda::EventAnalysis event;
    };
}

#endif
