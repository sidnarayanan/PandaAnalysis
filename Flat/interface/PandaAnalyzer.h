#ifndef PandaAnalyzer_h
#define PandaAnalyzer_h

#include "AnalyzerUtilities.h"
#include "GeneralTree.h"
#include "Selection.h"
#include "Module.h"
#include "Common.h"


namespace pa { 

    class PandaAnalyzer {
    public :

        PandaAnalyzer(Analysis* a, int debug_=0);
        ~PandaAnalyzer();
        void Reset();
        void Run();
        void Terminate();
        void AddPresel(Selection *s) { selections.push_back(s); }
        void AddGoodLumiRange(int run, int l0, int l1);

    private:
        bool PassGoodLumis(int run, int lumi);
        bool PassPresel(Selection::Stage stage);

        //////////////////////////////////////////////////////////////////////////////////////

        std::map<int,std::vector<LumiRange>> goodLumis;
        std::vector<Selection*> selections;

        // IO for the analyzer
        TFile *fOut{nullptr};     // output file is owned by PandaAnalyzer
        TTree *tOut{nullptr};
        GeneralTree gt; // essentially a wrapper around tOut
        TH1D *hDTotalMCWeight{nullptr};

        TFile *fIn{nullptr};
        TTree *tIn{nullptr};    // input tree to read
        panda::EventAnalysis event;

        std::vector<TString> wIDs;

        int DEBUG; //!< debug verbosity level
        Analysis *analysis{nullptr}; //!< configure what to run
        Registry registry;
        std::vector<BaseModule*> mods_reco, mods_gen, mods_all; 
        GlobalMod *gblmod;
        ConfigMod cfgmod; 

    };
}

#endif

