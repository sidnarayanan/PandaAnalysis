#include "../interface/Config.h"

using namespace pa;
using namespace std;

Utils::~Utils() 
{
  for (auto& f : fCorrs)  { if (f.get() != nullptr) f->Close(); }

  if (fMSDcorr.get())
    fMSDcorr->Close();
}

double Utils::getCorr(CorrectionType ct, double x, double y)
{
  if (f1Corrs[ct].get() != nullptr) {
    return f1Corrs[ct]->Eval(x);
  }

  if (h1Corrs[ct].get() != nullptr) {
    return h1Corrs[ct]->Eval(x);
  }

  if (h2Corrs[ct].get() != nullptr) {
    return h2Corrs[ct]->Eval(x,y);
  }
  return 0;
}

double Utils::getError(CorrectionType ct, double x, double y)
{
  if (f1Corrs[ct].get() != nullptr) {
    return 0;  
  }

  if (h1Corrs[ct].get() != nullptr) {
    return h1Corrs[ct]->Error(x);
  }

  if (h2Corrs[ct].get() != nullptr) {
    return h2Corrs[ct]->Error(x,y);
  }
  return 0;
}

void Utils::openCorr(CorrectionType ct, TString fpath, TString hname, int dim)
{
  fCorrs[ct].reset(TFile::Open(fpath));
  if (dim==1) 
    h1Corrs[ct].reset(new THCorr1(fCorrs[ct]->Get(hname)));
  else if (dim==2)
    h2Corrs[ct].reset(new THCorr2(fCorrs[ct]->Get(hname)));
  else 
    f1Corrs[ct].reset(new TF1Corr((TF1*)fCorrs[ct]->Get(hname))); 
}

