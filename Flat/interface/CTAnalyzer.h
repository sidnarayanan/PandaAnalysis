#ifndef CTAnalyzer_h
#define CTAnalyzer_h

#include "AnalyzerUtilities.h"
#include "Analyzer.h"
#include "JetGraphTree.h"
#include "CommonTopOp.h"
#include "PandaAnalysis/Flat/interface/Common.h"


namespace pa {

    class CTAnalyzer : public Analyzer<JetGraphTree> {
    public :
        CTAnalyzer(Analysis* a, int debug_=0);
        ~CTAnalyzer();
        void Reset();
        void Run();
        void Terminate();

    private:

        //////////////////////////////////////////////////////////////////////////////////////
        std::unique_ptr<CTOp> op{nullptr};
        CTEvent event;
    };
}

#endif
