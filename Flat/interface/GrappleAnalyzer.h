#ifndef GrappleAnalyzer_h
#define GrappleAnalyzer_h

#include "AnalyzerUtilities.h"
#include "Analyzer.h"
#include "EventTree.h"
#include "Selection.h"
#include "Operator.h"
#include "PandaAnalysis/Flat/interface/Common.h"

#include "CommonOps.h"

namespace pa {
    class GrappleOp : public EventOp {
    public:
        GrappleOp(panda::EventAnalysis& event_, Config& cfg_,
                  Utils& utils_, EventTree& gt_, int level_=0) :
            EventOp("grapple", event_, cfg_, utils_, gt_, level_) { }
        ~GrappleOp() { }
        virtual bool on() { return true; }

    protected:
        virtual void do_init(Registry& registry) {
            fOut = registry.access<TFile>("fOut");
            incrementAux(false);
        }
        virtual void do_execute();
        virtual void do_reset() {
            genJetInfo.reset();
        }
        virtual void do_terminate() {
            incrementAux(true);
        }
        void incrementAux(bool close = false);

    private:
        std::shared_ptr<TFile> fOut{nullptr};

        GenJetInfo genJetInfo;

        std::unique_ptr<TFile>fAux{nullptr};
        TTree* tAux{nullptr}; // don't ask why this is a bare pointer, I don't understand
                              // root memory management either.
        int auxCounter{0};
    };

    class GrappleAnalyzer : public Analyzer<EventTree> {
    public :
        GrappleAnalyzer(Analysis* a, int debug_=0);
        ~GrappleAnalyzer();
        void Reset();
        void Run();
        void Terminate();

    private:
        std::unique_ptr<GrappleOp> op{nullptr};
        Config cfg;
        Utils utils; 

        panda::EventAnalysis event;
    };
}

#endif
