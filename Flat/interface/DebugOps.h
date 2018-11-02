#include "Operator.h"

namespace pa {
    class PFDumpOp : public AnalysisOp {
        public:
            PFDumpOp(panda::EventAnalysis& event_,
                      Config& cfg_,
                      Utils& utils_,
                      GeneralTree& gt_,
                      int level_=0) : 
            AnalysisOp("pfdump", event_, cfg_, utils_, gt_, level_) { }
            virtual bool on() { return true; }
        protected:
            void do_execute() {
                logger.debug("PFDumpOp::do_execute", 
                             Form("event = %llu\n", gt.eventNumber));
                int i = 0;
                for (const auto& pf : event.pfCandidates) {
                    logger.debug(Form("PFDumpOp::do_execute %i", i),
                                 Form("  pt=%7.3f, eta=%6.3f, phi=%6.3f, pdgid=%6i", 
                                      pf.pt(), pf.eta(), pf.phi(), pf.pdgId()));
                    ++i; 
                }
            }
    };
}
