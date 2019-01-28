#ifndef JGAnalyzer_h
#define JGAnalyzer_h

#include "AnalyzerUtilities.h"
#include "Analyzer.h"
#include "JetGraphTree.h"
#include "Selection.h"
#include "Operator.h"
#include "PandaAnalysis/Flat/interface/Common.h"

#include "CommonOps.h"
#include "FatJetsOps.h"


namespace pa {

    class JGAnalyzer : public Analyzer<JetGraphTree> {
    public :
        JGAnalyzer(Analysis* a, int debug_=0);
        ~JGAnalyzer();
        void Reset();
        void Run();
        void Terminate();

    private:

        //////////////////////////////////////////////////////////////////////////////////////
        std::unique_ptr<BaseGenPOp<JetGraphTree>> gen{nullptr};
        std::unique_ptr<FatJetReclusterOp<JetGraphTree>> fj{nullptr};
        std::unique_ptr<JetGraphOp> op{nullptr};
        Config cfg;
        Utils utils; 

        panda::EventAnalysis event;
    };
}

#endif
