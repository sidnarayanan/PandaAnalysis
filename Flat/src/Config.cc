#include "../interface/Config.h"

using namespace pa;


Utils::~Utils() 
{
  for (auto* f : fCorrs)
    if (f)
      f->Close();

  for (auto* t : tCorrs)
    if (t)
      delete t; 

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
  delete rochesterCorrection;

  delete ecfcalc;
  delete grid;
}

double Utils::GetCorr(CorrectionType ct, double x, double y)
{
  TCorr *base = tCorr[ct];
  
  auto* fc = dynamic_cast<TFCorr*>(base);
  if (fc != nullptr) {
    return fc->Eval(x);
  }

  auto* h1c = dynamic_cast<THCorr1*>(base);
  if (h1c != nullptr) {
    return h1c->Eval(x);
  }

  auto* h2c = dynamic_cast<THCorr2*>(base);
  if (h2c != nullptr) {
    return h2c->Eval(x, y);
  }
}

double Utils::GetError(CorrectionType ct, double x, double y)
{
  TCorr *base = tCorr[ct];
  
  auto* h1c = dynamic_cast<THCorr1*>(base);
  if (h1c != nullptr) {
    return h1c->Error(x);
  }

  auto* h2c = dynamic_cast<THCorr2*>(base);
  if (h2c != nullptr) {
    return h2c->Error(x, y);
  }
}

void Utils::OpenCorrection(CorrectionType ct, TString fpath, TString hname, int dim)
{
  fCorrs[ct] = TFile::Open(fpath);
  if (dim==1) 
    h1Corrs[ct] = new THCorr1(fCorrs[ct]->Get(hname));
  else if (dim==2)
    h2Corrs[ct] = new THCorr2(fCorrs[ct]->Get(hname));
  else 
    f1Corrs[ct] = new TF1Corr((TF1*)fCorrs[ct]->Get(hname)); 
}