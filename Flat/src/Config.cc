#include "../interface/Config.h"

using namespace pa;
using namespace std;

Utils::~Utils() 
{
  for (auto* t : f1Corrs) { delete t; }
  for (auto* t : h1Corrs) { delete t; }
  for (auto* t : h2Corrs) { delete t; }
  for (auto* f : fCorrs)  { if (f != nullptr) f->Close(); }

  delete btag;
  delete eras; 

  delete activeArea;
  delete areaDef;

  if (fMSDcorr)
    fMSDcorr->Close();
}

double Utils::getCorr(CorrectionType ct, double x, double y)
{
  if (f1Corrs[ct] != nullptr) {
    return f1Corrs[ct]->Eval(x);
  }

  if (h1Corrs[ct] != nullptr) {
    return h1Corrs[ct]->Eval(x);
  }

  if (h2Corrs[ct] != nullptr) {
    return h2Corrs[ct]->Eval(x,y);
  }
  return 0;
}

double Utils::getError(CorrectionType ct, double x, double y)
{
  if (f1Corrs[ct] != nullptr) {
    return 0;  
  }

  if (h1Corrs[ct] != nullptr) {
    return h1Corrs[ct]->Error(x);
  }

  if (h2Corrs[ct] != nullptr) {
    return h2Corrs[ct]->Error(x,y);
  }
  return 0;
}

void Utils::openCorr(CorrectionType ct, TString fpath, TString hname, int dim)
{
  fCorrs[ct] = TFile::Open(fpath);
  if (dim==1) 
    h1Corrs[ct] = new THCorr1(fCorrs[ct]->Get(hname));
  else if (dim==2)
    h2Corrs[ct] = new THCorr2(fCorrs[ct]->Get(hname));
  else 
    f1Corrs[ct] = new TF1Corr((TF1*)fCorrs[ct]->Get(hname)); 
}

