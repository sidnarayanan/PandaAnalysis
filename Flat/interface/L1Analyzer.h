#ifndef L1Analyzer_h
#define L1Analyzer_h

#include "AnalyzerUtilities.h"
#include "Analyzer.h"
#include "L1Tree.h"
#include "L1Op.h"
#include "PandaAnalysis/Flat/interface/Common.h"


namespace pa {

    class L1Analyzer : public Analyzer<L1Tree> {
    public :
        L1Analyzer(Analysis* a, int debug_=0);
        ~L1Analyzer();
        void Reset();
        void Run();
        void Terminate();

    private:

        //////////////////////////////////////////////////////////////////////////////////////
        std::unique_ptr<L1Op> op{nullptr};
        L1Event event;
    };
}

#endif
