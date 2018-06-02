#ifndef KINFIT
#define KINFIT

#include "RooMinimizer.h"
#include "RooWorkspace.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "TLorentzVector.h"
#include "TString.h"
#include <vector>

namespace kinfit {


  class Fit {
  private:
    class Vector {
    public:
      Vector(int idx_) : 
        idx(idx_),
        pt(TString::Format("pt_%i", idx_).Data(),"",1),
        phi(TString::Format("phi_%i", idx_).Data(),"",1),
        scale(TString::Format("scale_%i", idx_).Data(),"",1, 0, 20),
        res(TString::Format("res_%i", idx_).Data(),"",1),
        px(TString::Format("px_%i", idx_).Data(),"",1),
        py(TString::Format("py_%i", idx_).Data(),"",1),
        pz(TString::Format("pz_%i", idx_).Data(),"",1),
        e(TString::Format("e_%i", idx_).Data(),"",1)
      {
      }

      void dump() { 
        printf("v%i -> (pt=%.3f, phi=%.3f, scale=%.3f, res=%.3f)\n", 
               idx, pt.getVal(), phi.getVal(), scale.getVal(), res.getVal());
      }

      int idx;
      TString name;
      RooRealVar pt, phi, scale, res, px, py, pz, e;
    };
  public:
    Fit(int N, float zMass = -1); 

    ~Fit() { delete fcn; delete px; 
             delete py; delete min; 
             delete constraint; delete mll; }

    void resetParticle(int idx) {
      auto& p = particles.at(idx);
      p.pt.setVal(0);
      p.phi.setVal(0);
      p.px.setVal(0);
      p.py.setVal(0);
      p.pz.setVal(0);
      p.e.setVal(0);
      p.scale.setVal(1); 
      p.res.setVal(0.00001);
    }
    void setParticle(int idx, const TLorentzVector& p4, float res, float scale=1) {
      auto& p = particles.at(idx);
      p.pt.setVal(p4.Pt());
      p.phi.setVal(p4.Phi());
      p.px.setVal(p4.Px());
      p.py.setVal(p4.Py());
      p.pz.setVal(p4.Pz());
      p.e.setVal(p4.E());
      p.scale.setVal(scale); 
      p.res.setVal(res);
    }

    float getScale(int idx) const {
      return particles.at(idx).scale.getVal();
    }

    void run(TString opts=""); 

    void setPrintLevel(int i) { min->setPrintLevel(i); }

  private:
    std::vector<Vector> particles;
    RooFormulaVar *fcn{nullptr}, *px{nullptr}, *py{nullptr}, *constraint{nullptr}, *mll{nullptr};
    RooMinimizer *min{nullptr};
  };

}

#endif
