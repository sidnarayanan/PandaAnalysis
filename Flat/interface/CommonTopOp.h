#include "Operator.h"

// this is treated separately from other Ops because
// it is not a real PandaTree-based Op

namespace pa {
  struct CTEvent {
    int is_signal_new;
    float E[200];
    float PX[200];
    float PY[200];
    float PZ[200];
  };

  class CTOp : public BaseOperator<JetGraphTree> {
  public:
    CTOp(JetGraphTree& gt_, CTEvent& event_) : BaseOperator<JetGraphTree>("CT", gt_), event(event_) { 
      jetDef.reset(new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.8));
    } 
    ~CTOp() { } 
    void reset() { gt.Reset(); } 
    void execute() {
      std::vector<fastjet::PseudoJet> particles;
      particles.reserve(100);
      for (int i = 0; i != 100; ++i) {
        particles.emplace_back(event.PX[i], event.PY[i], event.PZ[i], event.E[i]);
      }
      fastjet::ClusterSequence seq(particles, *jetDef);
      auto allJets = seq.inclusive_jets(0.);
      if (allJets.size() > 0) {

        auto& history = seq.history();
        auto&jets = seq.jets();
        std::map<const fastjet::PseudoJet*,int> jetIdx;

        gt.jetPdgId = event.is_signal_new; 

        for (auto& jet : jets) {
          if (jet.perp() > 0.01) {
            int idx = jetIdx.size();
            if (idx >= NNODE)
              break;
            jetIdx[&jet] = idx;
            gt.nodePt[idx] = jet.perp();
            gt.nodeEta[idx] = jet.eta();
            gt.nodePhi[idx] = jet.phi();
            gt.nodeE[idx] = jet.e();
          }
        }
        auto getIdx = [&](int i) { 
          if (i < 0)
            return -1;
          auto* addr = &(jets.at(i));
          if (jetIdx.find(addr) == jetIdx.end())
            return -1;
          int ret = jetIdx[addr];
          return (ret < NNODE) ? ret : -1;
        };

        for (auto& h : history) {
          int jetpIdx = h.jetp_index;
          jetpIdx = getIdx(jetpIdx);
          if (jetpIdx < 0)
            continue;

          gt.nodeIsRoot[jetpIdx] = (h.child < 0) ? 1 : 0;

          if (h.parent1 >= 0) {
            gt.nodeIsFinal[jetpIdx] = 0;
            auto& hp1 = history.at(h.parent1);
            int hp1Idx = getIdx(hp1.jetp_index);
            if (hp1Idx >= 0)
              gt.adj[jetpIdx][hp1Idx] = 1;
            auto& hp2 = history.at(h.parent2);
            int hp2Idx = getIdx(hp2.jetp_index);
            if (hp2Idx >= 0)
              gt.adj[jetpIdx][hp2Idx] = 1;
          } else {
            gt.nodeIsFinal[jetpIdx] = 1;
          }
        }

        gt.Fill();
      }
    } 
  protected:
    CTEvent& event; 
    std::unique_ptr<fastjet::JetDefinition> jetDef{nullptr};
  };
}
