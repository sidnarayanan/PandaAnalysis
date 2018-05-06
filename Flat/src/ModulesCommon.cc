#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>


using namespace panda;
using namespace std;

/***** miscellaneous helper functions *****/

bool PandaAnalyzer::PassPresel(Selection::Stage stage) 
{
  if (selections.size() == 0)
    return true;

  bool pass = false;
  for (auto* s : selections) {
    if (s->anded())
      continue;
    if (s->accept(stage)) {
      pass = true;
      break;
    }
  }

  for (auto* s : selections) {
    if (s->anded()) {
      pass = pass && s->accept(stage);
    }
  }

  return pass;
}


/***** common physics methods *****/

// Create a new auxillary file for reco info
// Responsible: S. Narayanan
/*
void PandaAnalyzer::IncrementAuxFile(bool close)
{
  if (fAux) {
    fAux->WriteTObject(tAux, "inputs", "Overwrite");
    fAux->Close();
  }
  if (close)
    return;

  TString path = TString::Format(auxFilePath.Data(),auxCounter++);
  fAux = TFile::Open(path.Data(), "RECREATE");
  if (DEBUG) PDebug("PandaAnalyzer::IncrementAuxFile", "Opening "+path);
  tAux = new TTree("inputs","inputs");
  
  pfInfo.resize(NMAXPF);
  for (int i = 0; i != NMAXPF; ++i) {
    pfInfo[i].resize(NPFPROPS);
  }
  tAux->Branch("kinematics",&pfInfo);
  
  svInfo.resize(NMAXSV);
  for (int i = 0; i != NMAXSV; ++i) {
    svInfo[i].resize(NSVPROPS);
  }
  tAux->Branch("svs",&svInfo);

  tAux->Branch("msd",&fjmsd,"msd/F");
  tAux->Branch("pt",&fjpt,"pt/F");
  tAux->Branch("rawpt",&fjrawpt,"rawpt/F");
  tAux->Branch("eta",&fjeta,"eta/F");
  tAux->Branch("phi",&fjphi,"phi/F");
  tAux->Branch("rho",&(gt->fjRho),"rho/f");
  tAux->Branch("rawrho",&(gt->fjRawRho),"rawrho/f");
  tAux->Branch("rho2",&(gt->fjRho2),"rho2/f");
  tAux->Branch("rawrho2",&(gt->fjRawRho2),"rawrho2/f");
  tAux->Branch("nPartons",&(gt->fjNPartons),"nPartons/I");
  tAux->Branch("nBPartons",&(gt->fjNBPartons),"nBPartons/I");
  tAux->Branch("nCPartons",&(gt->fjNCPartons),"nCPartons/I");
  tAux->Branch("partonM",&(gt->fjPartonM),"partonM/f");
  tAux->Branch("partonPt",&(gt->fjPartonPt),"partonPt/f");
  tAux->Branch("partonEta",&(gt->fjPartonEta),"partonEta/f");
  tAux->Branch("tau32",&(gt->fjTau32),"tau32/f");
  tAux->Branch("tau32SD",&(gt->fjTau32SD),"tau32SD/f");
  tAux->Branch("tau21",&(gt->fjTau21),"tau21/f");
  tAux->Branch("tau21SD",&(gt->fjTau21SD),"tau21SD/f");
  tAux->Branch("eventNumber",&(gt->eventNumber),"eventNumber/l");
  tAux->Branch("maxcsv",&(gt->fjMaxCSV),"maxcsv/f");
  tAux->Branch("mincsv",&(gt->fjMinCSV),"mincsv/f");
  tAux->Branch("doubleb",&(gt->fjDoubleCSV),"doubleb/f");

  gt->SetAuxTree(tAux);

  fOut->cd();

  if (tr)
    tr->TriggerEvent("increment aux file");
}

*/ 

