#ifndef PandaAnalyzer_h
#define PandaAnalyzer_h

#include "AnalyzerUtilities.h"
#include "GeneralTree.h"
#include "Selection.h"
#include "Module.h"
#include "PandaAnalysis/Flat/interface/Common.h"

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
        void AddPresel(Selection *s) { s->set_gt(&gt); selections.emplace_back(s); }
        void AddGoodLumiRange(int run, int l0, int l1);

        int firstEvent{0}, lastEvent{-1};
    private:
        bool PassGoodLumis(int run, int lumi);
        bool PassPresel(Selection::Stage stage);

        //////////////////////////////////////////////////////////////////////////////////////
        GeneralTree gt;
        int DEBUG; //!< debug verbosity level
        Analysis& analysis; //!< configure what to run
        Registry registry;

        std::vector<std::unique_ptr<AnalysisMod>> mods_all;
        GlobalMod *gblmod{nullptr};
        ContainerMod *preselmod{nullptr}, *postselmod{nullptr};
        ConfigMod cfgmod;

        std::map<int,std::vector<LumiRange>> goodLumis;
        std::vector<std::unique_ptr<Selection>> selections;

        // IO for the analyzer
        // output file is owned by PandaAnalyzer, but we use it elsewhere
        std::shared_ptr<TFile> fOut{nullptr};
        TTree *tOut{nullptr};

        std::unique_ptr<TFile> fIn{nullptr};
        TTree *tIn{nullptr};    // input tree to read
        panda::EventAnalysis event;

        std::shared_ptr<std::vector<TString>> wIDs;

    };
}

#endif
