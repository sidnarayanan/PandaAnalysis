#ifndef PandaAnalyzer_h
#define PandaAnalyzer_h

#include "AnalyzerUtilities.h"
#include "Analyzer.h"
#include "GeneralTree.h"
#include "Selection.h"
#include "Operator.h"
#include "PandaAnalysis/Flat/interface/Common.h"

#include "CommonOps.h"
#include "BTagOps.h"
#include "CommonOps.h"
#include "DeepOps.h"
#include "FatJetsOps.h"
#include "JetsOps.h"
#include "LepPhoOps.h"
#include "TheoryOps.h"
#include "HbbOps.h"
#include "DebugOps.h"
#include "CatOps.h"


namespace pa {

    class PandaAnalyzer : public Analyzer<GeneralTree> {
    public :
        PandaAnalyzer(Analysis* a, int debug_=0);
        ~PandaAnalyzer();
        void Reset();
        void Run();
        void Terminate();
        void AddPresel(Selection *s) { s->set_gt(&gt); selections.emplace_back(s); }
        void AddGoodLumiRange(int run, int l0, int l1);
    private:
        bool PassGoodLumis(int run, int lumi);
        bool PassPresel(Selection::Stage stage);

        //////////////////////////////////////////////////////////////////////////////////////
        std::vector<std::unique_ptr<AnalysisOp>> ops_all;
        GlobalOp *gblop{nullptr};
        ContainerOp *preselop{nullptr}, *postselop{nullptr};
        ConfigOp cfgop;

        std::map<int,std::vector<LumiRange>> goodLumis;
        std::vector<std::unique_ptr<Selection>> selections;

        panda::EventAnalysis event;

        std::shared_ptr<std::vector<TString>> wIDs;

    };
}

#endif
