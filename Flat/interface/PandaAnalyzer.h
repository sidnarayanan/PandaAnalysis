#ifndef PandaAnalyzer_h
#define PandaAnalyzer_h

#include "AnalyzerUtilities.h"
#include "GeneralTree.h"
#include "Selection.h"
#include "Module.h"
#include "Common.h"


namespace pa { 

    /////////////////////////////////////////////////////////////////////////////
    // PandaAnalyzer definition
    class PandaAnalyzer {
    public :

        PandaAnalyzer(int debug_=0);
        ~PandaAnalyzer();
        int Init(TTree *tree, TH1D *hweights, TTree *weightNames=0);
        void SetOutputFile(TString fOutName);
        void ResetBranches();
        void Run();
        void Terminate();
        void SetDataDir(const char *s);
        void AddPresel(Selection *s) { selections.push_back(s); }
        void AddGoodLumiRange(int run, int l0, int l1);

        // public configuration
        void SetAnalysis(Analysis *a) { analysis = a; }

    private:
        //////////////////////////////////////////////////////////////////////////////////////

        bool PassGoodLumis(int run, int lumi);
        bool PassPresel(Selection::Stage stage);

        //////////////////////////////////////////////////////////////////////////////////////

        int DEBUG = 0; //!< debug verbosity level
        Analysis *analysis{nullptr}; //!< configure what to run
        Registry registry;
        std::vector<BaseModule*> mods; 

        //////////////////////////////////////////////////////////////////////////////////////

        std::map<int,std::vector<LumiRange>> goodLumis;
        std::vector<Selection*> selections;

        // IO for the analyzer
        TString fOutPath;
        TFile *fOut{nullptr};     // output file is owned by PandaAnalyzer
        TTree *tOut{nullptr};
        GeneralTree gt; // essentially a wrapper around tOut
        TH1D *hDTotalMCWeight{nullptr};

        TTree *tIn{nullptr};    // input tree to read
        panda::EventAnalysis event;

        std::vector<TString> wIDs;
    };
}

#endif

