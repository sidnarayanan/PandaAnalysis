#include "PandaCore/Tools/interface/Common.h"
#include "AnalyzerUtilities.h"
#include "TString.h"
#include "GeneralTree.h"
#include <algorithm>

#ifndef PANDA_SELECTION
#define PANDA_SELECTION

#define __ACCPFUNC(x) [](const GeneralTree* gt) { return x; }

namespace pa {
  class Selection {
  public:
    enum Stage {
      sGen, sReco
    };

    Selection(Stage stage_, TString n=""): stage(stage_), name(n) { }
    virtual ~Selection() { }
    
    virtual void report() const final { 
      PDebug("Selection::" + name, Form("Accepted %i/%i events", nPassed, nTotal)); 
    }
    virtual void set_gt(const GeneralTree *gt_) final { gt = gt_; }
    // if called at a different stage, just return true
    virtual bool accept(Stage stage_) final { 
      if (stage_ != stage)
        return true;
      bool good = do_accept(); 
      ++nTotal; if (good) ++nPassed;
      return good;
    };
    virtual bool anded() const final { return is_anded; }
    virtual TString get_name() const final { return name; }
  protected:
    virtual bool do_accept() const = 0;
    Stage stage;
    int nTotal{0}, nPassed{0};
    const GeneralTree* gt{nullptr};
    TString name;
    bool is_anded{false};
  };

  // Lambda class that can be bound at runtime or subclassed (see below)
  class LambdaSel: public Selection {
  public:
    typedef std::function<bool(const GeneralTree*)> accept_func;
    LambdaSel(Stage stage_, TString n, accept_func f_, bool anded = false): 
      Selection(stage_,n), 
      f(f_) { is_anded = anded; }
    ~LambdaSel() { }

  protected:
    virtual bool do_accept() const final { return f(gt); }
  private:
    accept_func f;
  };

  class TriggerSel: public LambdaSel {
  public:
    TriggerSel():
      LambdaSel(Selection::sReco, "Trigger", 
                __ACCPFUNC((gt->isData==0) || (gt->trigger!=0)), 
                true) { }
  };

  class LowGenBosonPtSel: public LambdaSel {
  public:
    LowGenBosonPtSel():
      LambdaSel(Selection::sGen, "LowGenBosonPt", __ACCPFUNC(gt->trueGenBosonPt > 50 && gt->lheHT > 100)) { }
  };

  class GenBosonPtSel: public LambdaSel {
  public:
    GenBosonPtSel():
      LambdaSel(Selection::sGen, "GenBosonPt", __ACCPFUNC(gt->trueGenBosonPt > 100)) { }
  };

  class FatjetSel: public LambdaSel {
  public:
    FatjetSel():
      LambdaSel(Selection::sReco, "Fatjet", __ACCPFUNC(gt->fjPt[0] > 250)) { }
  };

  class Fatjet450Sel: public LambdaSel {
  public:
    Fatjet450Sel():
      LambdaSel(Selection::sReco, "Fatjet450", __ACCPFUNC(gt->fjPt[0] > 450)) { }
  };

  class GenFatJetSel: public LambdaSel {
  public:
    GenFatJetSel():
      LambdaSel(Selection::sGen, "GenFatJet", __ACCPFUNC(gt->genFatJetPt > 450)) { }
  };

  class GenFatJet200Sel: public LambdaSel {
  public:
    GenFatJet200Sel():
      LambdaSel(Selection::sGen, "GenFatJet200", __ACCPFUNC(gt->genFatJetPt > 200)) { }
  };


  class LeptonSel: public Selection {
  public:
    LeptonSel(): Selection(Selection::sReco, "lepton") { }

  protected:
    virtual bool do_accept() const;
  };


  class LeptonFakeSel: public Selection {
  public:
    LeptonFakeSel(): Selection(Selection::sReco, "lepton fake") { }

  protected:
    virtual bool do_accept() const; 
  };


  class RecoilSel: public Selection {
  public:
    RecoilSel(): Selection(Selection::sReco, "recoil") { }
    float threshold{200};
  protected:
    virtual bool do_accept() const;
    bool vary_jes{true};
  };
  typedef RecoilSel MonojetSel; // backwards compatibility


  class MonotopSel: public RecoilSel {
  public:
    MonotopSel(): RecoilSel() { vary_jes = false; name = "monotop"; }
  protected:
    virtual bool do_accept() const;
  };


  class MonohiggsSel: public RecoilSel {
  public:
    MonohiggsSel(): RecoilSel() { vary_jes = false; threshold = 175; name = "monohiggs"; }
  protected:
    virtual bool do_accept() const;
  };


  class VHbbSel: public Selection {
  public:
    VHbbSel(): Selection(Selection::sReco, "vhbb") { }
  protected:
    virtual bool do_accept() const;
  };
}

#endif
