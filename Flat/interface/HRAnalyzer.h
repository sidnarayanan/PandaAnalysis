#ifndef HRAnalyzer_h
#define HRAnalyzer_h

#include "AnalyzerUtilities.h"
#include "Analyzer.h"
#include "HeavyResTree.h"
#include "Selection.h"
#include "Module.h"
#include "PandaAnalysis/Flat/interface/Common.h"

#include "CommonMods.h"
#include "FatJetsMods.h"


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
        std::unique_ptr<HRGenPMod> gen{nullptr};
        std::unique_ptr<HRTagMod> mod{nullptr};
        Config cfg;
        Utils utils; 

        panda::EventAnalysis event;
    };
}

#endif
