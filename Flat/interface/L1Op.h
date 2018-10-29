#include "Operator.h"
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"

// this is treated separately from other Ops because
// it is not a real PandaTree-based Op

namespace pa {
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> LorentzVector;

  struct L1Event {
    Long64_t run;
    Long64_t lumi;
    Long64_t event;
    int bunchCrossing;
    int metFilter;

    LorentzVector *met_p4{nullptr};
    std::vector<LorentzVector> *jet_p4{nullptr};
    std::vector<float> *jet_neutralEmFrac{nullptr};
    std::vector<float> *jet_neutralHadFrac{nullptr};

    std::vector<int> *L1EG_bx{nullptr};
    std::vector<LorentzVector> *L1EG_p4{nullptr};
    std::vector<int> *L1EG_iso{nullptr};

    std::vector<int> *L1GtBx{nullptr};
  };

  class L1Op : public BaseOperator<L1Tree> {
  public:
    L1Op(L1Tree& gt_, L1Event& event_) : BaseOperator<L1Tree>("L1", gt_), event(event_) { } 
    ~L1Op() { } 
    void reset() { gt.Reset(); } 
    void execute(); 
  protected:
    L1Event& event; 
  };
}
