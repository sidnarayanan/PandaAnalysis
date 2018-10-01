#include "../interface/KinematicFit.h"

#define KINFITDEBUG 0


using namespace kinfit;


Fit::Fit(int N, float zMass) 
{
  particles.reserve(N);
  TString px_str{"0"}, py_str{"0"}, mll_str{"0"}, const_str{"0"};
  RooArgList args;
  static const int n_params = 8;
  for (int i = 0; i != N; ++i) {
    particles.emplace_back(i);
    auto& p = particles.back();
    args.add(static_cast<RooAbsArg&>(p.scale));
    args.add(static_cast<RooAbsArg&>(p.pt));
    args.add(static_cast<RooAbsArg&>(p.phi));
    args.add(static_cast<RooAbsArg&>(p.res));
    args.add(static_cast<RooAbsArg&>(p.e));
    args.add(static_cast<RooAbsArg&>(p.px));
    args.add(static_cast<RooAbsArg&>(p.py));
    args.add(static_cast<RooAbsArg&>(p.pz));
    px_str += TString::Format("+(@%i*@%i*cos(@%i))", n_params*i, n_params*i+1, n_params*i+2);
    py_str += TString::Format("+(@%i*@%i*sin(@%i))", n_params*i, n_params*i+1, n_params*i+2);
    const_str += TString::Format("+pow((@%i-1)/@%i,2)", n_params*i, n_params*i+3);
  }

  px = new RooFormulaVar("px", "px", px_str.Data(), args);
  py = new RooFormulaVar("py", "py", py_str.Data(), args);
  constraint = new RooFormulaVar("constraint", "constraint", const_str.Data(), args);
  if (N >= 2 && zMass > 0) {
    TString mll_str = TString::Format("pow((@%i*@%i)+(@%i*@%i),2)", 
                                      n_params*0+0, n_params*0+4, n_params*1+0, n_params*1+4);
    for (int param = 5; param != n_params; ++param) {
      mll_str += TString::Format("-pow((@%i*@%i)+(@%i*@%i),2)", 
                                 n_params*0+0, n_params*0+param, n_params*1+0, n_params*1+param);
    }
    mll_str += TString::Format(" - %f", std::pow(zMass,2));
    mll = new RooFormulaVar("mll", "mll", mll_str.Data(), args);
    // pT^2/(8 GeV)^2 + (mll^2)^2/Gamma_Z^4
    fcn = new RooFormulaVar("fcn", "fcn", "(pow(@0,2) + pow(@1,2))/64 + pow(@3,2)/39 + @2", 
                            RooArgList(*px, *py, *constraint, *mll));
  } else {
    fcn = new RooFormulaVar("fcn", "fcn", "(pow(@0,2) + pow(@1,2))/64 + @2", 
                            RooArgList(*px, *py, *constraint));
  }

  min = new RooMinimizer(dynamic_cast<RooAbsReal&>(*fcn));
//    min->setStrategy(2);
}



void Fit::run(TString opts) 
{
#if KINFITDEBUG
  printf("PX,PY,MASS,CONST = %f %f %f %f\n", px->getVal(), py->getVal(), mll->getVal(), constraint->getVal());
  for (auto& v : particles)
    v.dump();
#endif
  min->fit(opts.Data()); 
#if KINFITDEBUG
  printf("PX,PY,MASS,CONST = %f %f %f %f\n", px->getVal(), py->getVal(), mll->getVal(), constraint->getVal());
  for (auto& v : particles)
    v.dump();
#endif
}
