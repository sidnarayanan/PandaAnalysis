#ifndef L1Analyzer_h
#define L1Analyzer_h

#include "AnalyzerUtilities.h"
#include "L1Tree.h"
#include "L1Mod.h"
#include "PandaAnalysis/Flat/interface/Common.h"


namespace pa {

    class L1Analyzer {
    public :
        L1Analyzer(Analysis* a, int debug_=0);
        ~L1Analyzer();
        void Reset();
        void Run();
        void Terminate();

        int firstEvent{0}, lastEvent{-1};
    private:

        //////////////////////////////////////////////////////////////////////////////////////
        L1Tree gt;
        int DEBUG; //!< debug verbosity level
        Analysis& analysis; 
        Registry registry;

        std::unique_ptr<L1Mod> mod{nullptr};

        // IO for the analyzer
        // output file is owned by L1Analyzer, but we use it elsewhere
        std::shared_ptr<TFile> fOut{nullptr};
        TTree *tOut{nullptr};

        std::unique_ptr<TFile> fIn{nullptr};
        TTree *tIn{nullptr};    // input tree to read
        L1Event event;
    };
}

#endif
