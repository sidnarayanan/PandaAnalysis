#ifndef HRAnalyzer_h
#define HRAnalyzer_h

#include "AnalyzerUtilities.h"
#include "HeavyResTree.h"
#include "Selection.h"
#include "Module.h"
#include "PandaAnalysis/Flat/interface/Common.h"

#include "CommonMods.h"
#include "FatJetsMods.h"


namespace pa {

    class HRAnalyzer {
    public :
        HRAnalyzer(Analysis* a, int debug_=0);
        ~HRAnalyzer();
        void Reset();
        void Run();
        void Terminate();

        int firstEvent{0}, lastEvent{-1};
    private:

        //////////////////////////////////////////////////////////////////////////////////////
        HeavyResTree gt;
        int DEBUG; //!< debug verbosity level
        Analysis& analysis; //!< configure what to run
        Registry registry;

        std::unique_ptr<HRGenPMod> gen{nullptr};
        std::unique_ptr<HRTagMod> mod{nullptr};
        Config cfg;
        Utils utils; 

        // IO for the analyzer
        // output file is owned by HRAnalyzer, but we use it elsewhere
        std::shared_ptr<TFile> fOut{nullptr};
        TTree *tOut{nullptr};

        std::unique_ptr<TFile> fIn{nullptr};
        TTree *tIn{nullptr};    // input tree to read
        panda::EventAnalysis event;
    };
}

#endif
