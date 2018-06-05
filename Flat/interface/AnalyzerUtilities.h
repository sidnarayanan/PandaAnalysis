#ifndef ANALYZERUTILS_H
#define ANALYZERUTILS_H

// fastjet
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"

// root, stl
#include "TRotation.h"
#include "TSystem.h"
#include <set>
#include <stdexcept>

#include "PandaAnalysis/Flat/interface/Common.h"

namespace std{
  template <typename T>
    inline shared_ptr<vector<T>> v_make_shared(int n = 0) { return make_shared<vector<T>>(n); }
};

namespace pa { 
  void downloadData(TString url, TString outpath, bool force=true, TString opts="-nv"); 

  // semi-temporary measure to deal with v009 gen duplication issue
  inline const panda::GenParticle* pToGPtr(const panda::Particle* p)
  {
    auto* g = dynamic_cast<const panda::GenParticle*>(p);
    if (g == nullptr) {
      logger.error("pToGPtr",Form("Trying to cast a non-GenParticle at %p", p));
      throw std::runtime_error("");
    }
    return g;
  }

  inline const panda::GenParticle& pToGRef(const panda::Particle* p)
  {
    // I think this is what I want - the underlying data type really is a GenParticle
    return *pToGPtr(p);
  }

  ////////////////////////////////////////////////////////////////////////////////////

  class JetTree {
    public:
      JetTree(fastjet::PseudoJet& root): _root(root) { }
      ~JetTree() { }
      std::vector<int> GetTerminals() {
        std::vector<int> v;
        _root.GetTerminals(v);
        return v;
      }

      struct compare : public std::unary_function<fastjet::PseudoJet, bool> {
        explicit compare(int x_) : _x(x_) {}
        bool operator() (const fastjet::PseudoJet& y_) { return _x == y_.user_index(); }
        int _x;
      };

    private:
      class Node {
        public:
          Node(fastjet::PseudoJet& pj_);
          ~Node() { delete l; delete r; }
          void GetTerminals(std::vector<int>&);
          fastjet::PseudoJet _pj;
          Node *l=0, *r=0;   
      };

      Node _root;
  };

  ////////////////////////////////////////////////////////////////////////////////////

  class ParticleGridder {
    public:
      ParticleGridder(unsigned etaN, unsigned phiN, float etaMax=5);
      ~ParticleGridder() { clear(); }
      void clear();
      void add(const panda::Particle& p);
      std::vector<TLorentzVector>& get();
      bool _on = true;
      bool _etaphi = true;
    private:
      class Bin {
      public:
        Bin(float lo, float hi, int N):
          _lo(lo),
          _scale(2. * N / (hi - lo)),
          _invscale(1. / _scale)
        {}
        int find(float x) { return static_cast<int>( (x - _lo) * _scale ); }
        float center(int i) { return (_invscale * (i + 0.5)) + _lo; }
        float left(int i) { return (_invscale * i) + _lo; }
      private:
        float _lo, _scale, _invscale;
      }; 

      float _etaMax, _phiMax;
      Bin _etaBin, _phiBin;
      std::map<std::pair<int,int>, std::vector<TLorentzVector*>> _collections;
      std::vector<TLorentzVector> _particles;
      std::vector<TLorentzVector> _gridded; 
  };

  ////////////////////////////////////////////////////////////////////////////////////

  class JetRotation {
    public:
      JetRotation(float x1, float y1, float z1,
                  float x2, float y2, float z2);
      void Rotate(float& x, float& y, float& z);
    private:
      TRotation r_toz;
      TRotation r_inxy;
  };

  ////////////////////////////////////////////////////////////////////////////////////

  typedef std::vector<fastjet::PseudoJet> VPseudoJet;
  VPseudoJet convertPFCands(const std::vector<const panda::PFCand*> &incoll, bool puppi, double minPt=0.001);
  VPseudoJet convertPFCands(const panda::RefVector<panda::PFCand> &incoll, bool puppi, double minPt=0.001);
  VPseudoJet convertPFCands(const panda::PFCandCollection &incoll, bool puppi, double minPt=0.001);

  ////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////

  class LumiRange {
  public:
      LumiRange(int l0_,int l1_):
          l0(l0_),
          l1(l1_)
       { }
      ~LumiRange() {}
      bool Contains(int l) {
          return l0<=l && l<=l1;
      }
  private:
      int l0, l1;
  };

  ////////////////////////////////////////////////////////////////////////////////////

  class TriggerHandler {  
  public:
    TriggerHandler() {};
    ~TriggerHandler() {};
    void addTriggers(std::vector<TString> paths_) { 
      for (auto &path : paths_) {
        paths.push_back(path); 
        indices.push_back(-1); 
      }
    }
    void registryTrigger(unsigned my_idx, int panda_idx) { indices[my_idx] = panda_idx; }
    std::vector<int> indices;
    std::vector<TString> paths;
  };

  ////////////////////////////////////////////////////////////////////////////////////


  class TF1Corr {
  public:
    TF1Corr(TF1 *f_) {
      h = f_;
    }
    ~TF1Corr() {} 
    double Eval(double x) {
      return h->Eval(x);
    }

    TF1 *GetFunc() { return h; }
  private:
    TF1* h=nullptr;
  };

  template <typename T>
  class THCorr {
    public:
      // wrapper around TH[12] to do corrections
      THCorr(TObject *h_) {
        this->h = new T();
        this->h->SetDirectory(0);
        h_->Copy(*(this->h));  // easiest way to cast e.g. TH1F->TH1D
        dim = this->h->GetDimension();
        TAxis *thurn = this->h->GetXaxis(); 
        lo1 = thurn->GetBinCenter(1);
        hi1 = thurn->GetBinCenter(thurn->GetNbins());
        if (dim>1) {
          TAxis *taxis = this->h->GetYaxis();
          lo2 = taxis->GetBinCenter(1);
          hi2 = taxis->GetBinCenter(taxis->GetNbins());
        }
      }
      ~THCorr() { /* delete this->h;*/ } // ROOT is doing something to this memory that I don't understand
                                         // code as-is is causing a small memory leak 
      double Eval(double x) {
        if (dim!=1) {
          logger.error("THCorr1::Eval",
            TString::Format("Trying to access a non-1D histogram (%s)!",this->h->GetName()));
          return -1;
        }
        return getVal(this->h,bound(x,lo1,hi1));
      }

      double Eval(double x, double y) {
        if (dim!=2) {
          logger.error("THCorr1::Error",
           TString::Format("Trying to access a non-2D histogram (%s)!",this->h->GetName()));
          return -1;
        }
        return getVal(this->h,bound(x,lo1,hi1),bound(y,lo2,hi2));
      }

      double Error(double x) {
        if (dim!=1) {
          logger.error("THCorr1::Error",
            TString::Format("Trying to access a non-1D histogram (%s)!",this->h->GetName()));
          return -1;
        }
        return getError(this->h,bound(x,lo1,hi1));
      }

      double Error(double x, double y) {
        if (dim!=2) {
          logger.error("THCorr1::Eval",
           TString::Format("Trying to access a non-2D histogram (%s)!",this->h->GetName()));
          return -1;
        }
        return getError(this->h,bound(x,lo1,hi1),bound(y,lo2,hi2));
      }

      T *GetHist() { return this->h; }

    private:
      int dim;
      double lo1, lo2, hi1, hi2;
      T* h{nullptr};
  };

  typedef THCorr<TH1D> THCorr1;
  typedef THCorr<TH2D> THCorr2;

  ////////////////////////////////////////////////////////////////////////////////////

  bool ElectronIP(double eta, double dxy, double dz); 
  bool MuonIP(double dxy, double dz); 

  template <typename T>
  bool isMatched(const std::vector<T*>*objects,
                 double deltaR2, double eta, double phi) {
    for (auto *x : *objects) {
      if (x->pt()>0) {
        if ( DeltaR2(x->eta(),x->phi(),eta,phi) < deltaR2 )
          return true;
      }
    }
    return false;
  }

////////////////////////////////////////////////////////////////////////////////////
  struct GenJetInfo {
  public:
    float pt=-1, eta=-1, phi=-1, m=-1;
    float msd=-1;
    float tau3=-1, tau2=-1, tau1=-1;
    float tau3sd=-1, tau2sd=-1, tau1sd=-1;
    int nprongs=-1;
    float partonpt=-1, partonm=-1;
    std::vector<std::vector<float>> particles;
    std::vector<std::vector<std::vector<float>>> ecfs; // uh
    void reset() {
    pt=-1; eta=-1; phi=-1; m=-1;
    msd=-1;
    tau3=-1; tau2=-1; tau1=-1;
    tau3sd=-1; tau2sd=-1; tau1sd=-1;
    nprongs=-1;
    partonpt=-1; partonm=-1;
    for (auto& v : particles) {
      std::fill(v.begin(), v.end(), 0);
    }
    for (auto& v : ecfs) {
      for (auto& vv : v) {
      std::fill(vv.begin(), vv.end(), -1);
      }
    }
    }
  };

  struct JetHistory {
    int user_idx;
    int child_idx;
  };

  struct JetWrapper {
    float pt; 
    int flavor{0};
    float genpt{0};
    float maxpt{0}; 
    bool iso{false};
    int user_idx{-1};
    float breg{-1}, bregwidth{-1}; 
    const panda::Jet* base;
    const JetWrapper* nominal{nullptr}; 

    JetWrapper(float pt_, const panda::Jet& j): pt(pt_), base(&j) { }

    const panda::Jet& get_base() const { return *base; }
    float scale() const { return pt / base->pt(); }
    void p4(TLorentzVector& v) const {
      v.SetPtEtaPhiM(pt, base->eta(), base->phi(), base->m()); 
    }
    TLorentzVector p4() const { 
      TLorentzVector v; 
      p4(v);
      return v;
    }

  };

  struct JESHandler {
    std::vector<JetWrapper> all; // all jets 
    std::vector<JetWrapper*> cleaned;
    std::vector<JetWrapper*> cleaned_sorted; // pt-sorted
    std::vector<JetWrapper*> iso;   // cleaned that do not overlap with fj
    std::vector<JetWrapper*> central; // cleaned that are central
    std::vector<JetWrapper*> bcand;   // all that are b candidates

    TLorentzVector vpfMET, vpuppiMET;
    TVector2 vpfMETNoMu;
    TLorentzVector vpfUW, vpfUZ, vpfUA;
    TLorentzVector vpuppiUW, vpuppiUZ, vpuppiUA;
    int shift_idx{-1};

    void sort() { 
      for (auto* jw : cleaned)
        cleaned_sorted.push_back(jw);
      std::sort(cleaned_sorted.begin(), cleaned_sorted.end(), 
                [](const JetWrapper* x, const JetWrapper* y) { return x->pt > y->pt; });
    }
    void clear() { 
        all.clear(); cleaned.clear(); iso.clear();
        cleaned_sorted.clear();
        central.clear(); bcand.clear(); 
        vpfMETNoMu.SetMagPhi(0,0);
        for (TLorentzVector* v_ : {&vpfMET, &vpuppiMET, 
                     &vpfUW, &vpfUZ, &vpfUA,
                     &vpuppiUW, &vpuppiUZ, &vpuppiUA})
          v_->SetPtEtaPhiM(0,0,0,0);
      }
    void reserve(int N) { all.reserve(N); cleaned.reserve(N); iso.reserve(N);
              cleaned_sorted.reserve(N);
              central.reserve(N); bcand.reserve(N); }
  };
}

#endif
