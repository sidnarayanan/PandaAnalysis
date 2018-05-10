#ifndef PandaAnalyzer_h
#define PandaAnalyzer_h

#include "AnalyzerUtilities.h"
#include "GeneralTree.h"
#include "Selection.h"
#include "Module.h"
#include "Common.h"

#include "CommonMods.h"
#include "BTagMods.h"
#include "CommonMods.h"
#include "DeepMods.h"
#include "FatJetsMods.h"
#include "JetsMods.h"
#include "LepPhoMods.h"
#include "TheoryMods.h"
#include "HbbMods.h"


namespace pa {

    class PandaAnalyzer {
    public :
        PandaAnalyzer(Analysis* a, int debug_=0);
        ~PandaAnalyzer();
        void Reset();
        void Run();
        void Terminate();
        void AddPresel(Selection *s) { s->set_gt(&gt); selections.push_back(s); }
        void AddGoodLumiRange(int run, int l0, int l1);

        int firstEvent{0}, lastEvent{-1};
    private:
        bool PassGoodLumis(int run, int lumi);
        bool PassPresel(Selection::Stage stage);

        //////////////////////////////////////////////////////////////////////////////////////

        std::map<int,std::vector<LumiRange>> goodLumis;
        std::vector<Selection*> selections;

        // IO for the analyzer
        TFile *fOut{nullptr};     // output file is owned by PandaAnalyzer
        TTree *tOut{nullptr};
        GeneralTree gt;

        TFile *fIn{nullptr};
        TTree *tIn{nullptr};    // input tree to read
        panda::EventAnalysis event;

        std::vector<TString> wIDs;

        int DEBUG; //!< debug verbosity level
        Analysis& analysis; //!< configure what to run
        Registry registry;
        std::vector<AnalysisMod*> mods_all;
        GlobalMod *gblmod{nullptr};
        ContainerMod *preselmod{nullptr}, *postselmod{nullptr};
        ConfigMod cfgmod;

    };
}

#endif
