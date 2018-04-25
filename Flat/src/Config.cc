#include "../interface/Config.h"

using namespace pa;
using namespace std;

Utils::~Utils() 
{
  for (auto* t : f1Corrs) { if (t != nullptr) delete t; }
  for (auto* t : h1Corrs) { if (t != nullptr) delete t; }
  for (auto* t : h2Corrs) { if (t != nullptr) delete t; }
  for (auto* f : fCorrs)  { if (f != nullptr) f->Close(); }

  delete btagCalib;
  delete sj_btagCalib;
  for (auto* reader : btagReaders)
    delete reader;

  for (auto& iter : ak8UncReaders)
    delete iter.second;

  delete ak8JERReader;

  for (auto& iter : ak4UncReaders)
    delete iter.second;

  for (auto& iter : ak4ScaleReader)
    delete iter.second;

  delete ak4JERReader;

  delete activeArea
  delete areaDef;
  delete jetDef;
  delete jetDefKt;
  delete jetDefGen;
  delete softDrop;

  delete hDTotalMCWeight;
  
  delete bjetregReader;
  delete bjetreg_vars;
  delete rochesterCorrection;

  delete ecfcalc;
  delete grid;

  delete rochesterCorrection;
  delete csvReweighter;
  delete cmvaReweighter;
  if (fMSDCorr)
    fMSDCorr->Close();
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

}

double Utils::getError(CorrectionType ct, double x, double y)
{
  if (f1Corrs[ct] != nullptr) {
    return f1Corrs[ct]->Error(x);
  }

  if (h1Corrs[ct] != nullptr) {
    return h1Corrs[ct]->Error(x);
  }

  if (h2Corrs[ct] != nullptr) {
    return h2Corrs[ct]->Error(x,y);
  }
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

