#include "../interface/PandaLeptonicAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

#include "TSystem.h"
#include "TRandom2.h"

#include "PandaAnalysis/Flat/interface/LeptonicFactors.h"

using namespace panda;
using namespace std;

PandaLeptonicAnalyzer::PandaLeptonicAnalyzer(int debug_/*=0*/) {
  DEBUG = debug_;

  if (DEBUG) PDebug("PandaLeptonicAnalyzer::PandaLeptonicAnalyzer","Calling constructor");
  gt = new GeneralLeptonicTree();
  if (DEBUG) PDebug("PandaLeptonicAnalyzer::PandaLeptonicAnalyzer","Built GeneralLeptonicTree");
  flags["puppi"]     = false;
  flags["applyJSON"] = true;
  flags["genOnly"]   = false;
  flags["lepton"]    = false;
  if (DEBUG) PDebug("PandaLeptonicAnalyzer::PandaLeptonicAnalyzer","Called constructor");
}


PandaLeptonicAnalyzer::~PandaLeptonicAnalyzer() {
  if (DEBUG) PDebug("PandaLeptonicAnalyzer::~PandaLeptonicAnalyzer","Calling destructor");
}


void PandaLeptonicAnalyzer::ResetBranches() {
  genObjects.clear();
  matchPhos.clear();
  matchEles.clear();
  matchLeps.clear();
  gt->Reset();
  if (DEBUG) PDebug("PandaLeptonicAnalyzer::ResetBranches","Reset");
}


void PandaLeptonicAnalyzer::SetOutputFile(TString fOutName) {
  fOut = new TFile(fOutName,"RECREATE");
  tOut = new TTree("events","events");

  fOut->WriteTObject(hDTotalMCWeight);

  // fill the signal weights
  for (auto& id : wIDs) 
    gt->signal_weights[id] = 1;

  // Build the input tree here 
  gt->WriteTree(tOut);

  if (DEBUG) PDebug("PandaLeptonicAnalyzer::SetOutputFile","Created output in "+fOutName);
}


// phi*
double PandaLeptonicAnalyzer::phi_star_eta(TLorentzVector lep1, TLorentzVector lep2, int pdgId1) {
double theta_star_eta = 0;
if(pdgId1 > 0){ // pdgId > 0 == q < 0
  theta_star_eta = TMath::ACos(TMath::TanH((lep1.Eta()-lep2.Eta())/2.0));
} else {
  theta_star_eta = TMath::ACos(TMath::TanH((lep2.Eta()-lep1.Eta())/2.0));
}

double dphi = TMath::Abs(lep1.DeltaPhi(lep2));

return TMath::Tan((TMath::Pi()-dphi)/2.0) * TMath::Sin(theta_star_eta);
}

int PandaLeptonicAnalyzer::Init(TTree *t, TH1D *hweights, TTree *weightNames)
{
  if (DEBUG) PDebug("PandaLeptonicAnalyzer::Init","Starting initialization");
  if (!t || !hweights) {
    PError("PandaLeptonicAnalyzer::Init","Malformed input!");
    return 0;
  }
  tIn = t;

  event.setStatus(*t, {"!*"}); // turn everything off first

  TString jetname = (flags["puppi"]) ? "puppi" : "chs";
  panda::utils::BranchList readlist({"runNumber", "lumiNumber", "eventNumber", "rho", 
                                     "isData", "npv", "npvTrue", "weight", "chsAK4Jets", 
                                     "electrons", "muons", "taus", "photons", 
                                     "pfMet", "caloMet", "puppiMet", "rawMet", "trkMet",
                                     "recoil","metFilters","genMet","superClusters", "vertices", "triggerObjects"});
  readlist.setVerbosity(0);

  readlist.push_back("triggers");

  if (isData) {
  } else {
   readlist.push_back("genParticles");
   readlist.push_back("genReweight");
   readlist.push_back("partons");
   readlist.push_back("ak4GenJets");
  }

  event.setAddress(*t, readlist); // pass the readlist so only the relevant branches are turned on
  if (DEBUG) PDebug("PandaLeptonicAnalyzer::Init","Set addresses");

  hDTotalMCWeight = new TH1F("hDTotalMCWeight","hDTotalMCWeight",1,0,2);
  hDTotalMCWeight->SetBinContent(1,hweights->GetBinContent(1));

  Float_t xbinsPt[nBinPt+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,400,500,800,1500};
  Float_t xbinsRap[nBinRap+1] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4};
  Float_t xbinsPhiStar[nBinPhiStar+1] = {1e-3, 2e-3, 3e-3, 4e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3,
                                         1e-2, 2e-2, 3e-2, 4e-2, 5e-2, 6e-2, 7e-2, 8e-2, 9e-2,
                                         1e-1, 2e-1, 3e-1, 4e-1, 5e-1, 6e-1, 7e-1, 8e-1, 9e-1,
                                            1,    3,    5,    7,   10,   20,   30,   50};
  Float_t xbinsPtRap0[nBinPtRap0+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,400,1500};
  Float_t xbinsPtRap1[nBinPtRap1+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,400,1500};
  Float_t xbinsPtRap2[nBinPtRap2+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,400,1500};
  Float_t xbinsPtRap3[nBinPtRap3+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,400,1500};
  Float_t xbinsPtRap4[nBinPtRap4+1] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,20,22,25,28,32,37,43,52,65,85,120,160,190,220,250,300,400,1500};

  Float_t xbinsWWSSMLL[nBinWWSS+1];
    xbinsWWSSMLL[ 0] =   0;      xbinsWWSSMLL[ 1] = 100;      xbinsWWSSMLL[ 2] = 200;      xbinsWWSSMLL[ 3] = 300;      xbinsWWSSMLL[ 4] = 400;
    xbinsWWSSMLL[ 5] = 600;

  Float_t xbinsWWMLL[nBinWW+1];
    xbinsWWMLL[ 0] =  20;      xbinsWWMLL[ 1] =  40;      xbinsWWMLL[ 2] =  55;      xbinsWWMLL[ 3] =  70;      xbinsWWMLL[ 4] =  85;
    xbinsWWMLL[ 5] = 100;      xbinsWWMLL[ 6] = 120;      xbinsWWMLL[ 7] = 150;      xbinsWWMLL[ 8] = 200;      xbinsWWMLL[ 9] = 800;

  Float_t xbinsWWPTL1[nBinWW+1];
    xbinsWWPTL1[ 0] =  25;      xbinsWWPTL1[ 1] =  35;      xbinsWWPTL1[ 2] =  40;      xbinsWWPTL1[ 3] =  45;      xbinsWWPTL1[ 4] =  50;
    xbinsWWPTL1[ 5] =  60;      xbinsWWPTL1[ 6] =  70;      xbinsWWPTL1[ 7] =  80;      xbinsWWPTL1[ 8] = 110;      xbinsWWPTL1[ 9] = 150;

  Float_t xbinsWWPTL2[nBinWW+1];
    xbinsWWPTL2[ 0] =  20;      xbinsWWPTL2[ 1] =  25;      xbinsWWPTL2[ 2] =  30;      xbinsWWPTL2[ 3] =  35;      xbinsWWPTL2[ 4] =  40;
    xbinsWWPTL2[ 5] =  45;      xbinsWWPTL2[ 6] =  50;      xbinsWWPTL2[ 7] =  75;      xbinsWWPTL2[ 8] = 100;      xbinsWWPTL2[ 9] = 150;

  Float_t xbinsWWDPHILL[nBinWW+1];
    xbinsWWDPHILL[ 0] =   0*TMath::Pi()/180.;      xbinsWWDPHILL[ 1] =  20*TMath::Pi()/180.;      xbinsWWDPHILL[ 2] =  40*TMath::Pi()/180.;      xbinsWWDPHILL[ 3] =  60*TMath::Pi()/180.;      xbinsWWDPHILL[ 4] =  80*TMath::Pi()/180.;
    xbinsWWDPHILL[ 5] = 100*TMath::Pi()/180.;      xbinsWWDPHILL[ 6] = 120*TMath::Pi()/180.;      xbinsWWDPHILL[ 7] = 140*TMath::Pi()/180.;      xbinsWWDPHILL[ 8] = 160*TMath::Pi()/180.;      xbinsWWDPHILL[ 9] = 180*TMath::Pi()/180.;

  Float_t xbinsWWN0JET[nBinWWN0JET+1];
    xbinsWWN0JET[ 0] =-0.5;      xbinsWWN0JET[ 1] = 0.5;      xbinsWWN0JET[ 2] = 1.5;      xbinsWWN0JET[ 3] = 2.5;

  hDDilPtMM = new TH1D("hDDilPtMM", "hDDilPtMM", nBinPt, xbinsPt);
  hDDilPtEE = new TH1D("hDDilPtEE", "hDDilPtEE", nBinPt, xbinsPt);
  hDDilRapMM = new TH1D("hDDilRapMM", "hDDilRapMM", nBinRap, xbinsRap);
  hDDilRapEE = new TH1D("hDDilRapEE", "hDDilRapEE", nBinRap, xbinsRap);
  hDDilPhiStarMM = new TH1D("hDDilPhiStarMM", "hDDilPhiStarMM", nBinPhiStar, xbinsPhiStar);
  hDDilPhiStarEE = new TH1D("hDDilPhiStarEE", "hDDilPhiStarEE", nBinPhiStar, xbinsPhiStar);
  hDDilPtRap0MM = new TH1D("hDDilPtRap0MM", "hDDilPtRap0MM", nBinPtRap0, xbinsPtRap0);
  hDDilPtRap0EE = new TH1D("hDDilPtRap0EE", "hDDilPtRap0EE", nBinPtRap0, xbinsPtRap0);
  hDDilPtRap1MM = new TH1D("hDDilPtRap1MM", "hDDilPtRap1MM", nBinPtRap1, xbinsPtRap1);
  hDDilPtRap1EE = new TH1D("hDDilPtRap1EE", "hDDilPtRap1EE", nBinPtRap1, xbinsPtRap1);
  hDDilPtRap2MM = new TH1D("hDDilPtRap2MM", "hDDilPtRap2MM", nBinPtRap2, xbinsPtRap2);
  hDDilPtRap2EE = new TH1D("hDDilPtRap2EE", "hDDilPtRap2EE", nBinPtRap2, xbinsPtRap2);
  hDDilPtRap3MM = new TH1D("hDDilPtRap3MM", "hDDilPtRap3MM", nBinPtRap3, xbinsPtRap3);
  hDDilPtRap3EE = new TH1D("hDDilPtRap3EE", "hDDilPtRap3EE", nBinPtRap3, xbinsPtRap3);
  hDDilPtRap4MM = new TH1D("hDDilPtRap4MM", "hDDilPtRap4MM", nBinPtRap4, xbinsPtRap4);
  hDDilPtRap4EE = new TH1D("hDDilPtRap4EE", "hDDilPtRap4EE", nBinPtRap4, xbinsPtRap4);
  hDWWSSMLL = new TH1D("hDWWSSMLL", "hDWWSSMLL", nBinWWSS, xbinsWWSSMLL);
  hDWWMLL = new TH1D("hDWWMLL", "hDWWMLL", nBinWW, xbinsWWMLL);
  hDWWPTL1 = new TH1D("hDWWPTL1", "hDWWPTL1", nBinWW, xbinsWWPTL1);
  hDWWPTL2 = new TH1D("hDWWPTL2", "hDWWPTL2", nBinWW, xbinsWWPTL2);
  hDWWDPHILL = new TH1D("hDWWDPHILL", "hDWWDPHILL", nBinWW, xbinsWWDPHILL);
  hDWWN0JET = new TH1D("hDWWN0JET", "hDWWN0JET", nBinWWN0JET, xbinsWWN0JET);

  hDDilPtMM_PDF = new TH1D("hDDilPtMM_PDF", "hDDilPtMM_PDF", nBinPt, xbinsPt);
  hDDilPtEE_PDF = new TH1D("hDDilPtEE_PDF", "hDDilPtEE_PDF", nBinPt, xbinsPt);
  hDDilRapMM_PDF = new TH1D("hDDilRapMM_PDF", "hDDilRapMM_PDF", nBinRap, xbinsRap);
  hDDilRapEE_PDF = new TH1D("hDDilRapEE_PDF", "hDDilRapEE_PDF", nBinRap, xbinsRap);
  hDDilPhiStarMM_PDF = new TH1D("hDDilPhiStarMM_PDF", "hDDilPhiStarMM_PDF", nBinPhiStar, xbinsPhiStar);
  hDDilPhiStarEE_PDF = new TH1D("hDDilPhiStarEE_PDF", "hDDilPhiStarEE_PDF", nBinPhiStar, xbinsPhiStar);
  hDDilPtRap0MM_PDF = new TH1D("hDDilPtRap0MM_PDF", "hDDilPtRap0MM_PDF", nBinPtRap0, xbinsPtRap0);
  hDDilPtRap0EE_PDF = new TH1D("hDDilPtRap0EE_PDF", "hDDilPtRap0EE_PDF", nBinPtRap0, xbinsPtRap0);
  hDDilPtRap1MM_PDF = new TH1D("hDDilPtRap1MM_PDF", "hDDilPtRap1MM_PDF", nBinPtRap1, xbinsPtRap1);
  hDDilPtRap1EE_PDF = new TH1D("hDDilPtRap1EE_PDF", "hDDilPtRap1EE_PDF", nBinPtRap1, xbinsPtRap1);
  hDDilPtRap2MM_PDF = new TH1D("hDDilPtRap2MM_PDF", "hDDilPtRap2MM_PDF", nBinPtRap2, xbinsPtRap2);
  hDDilPtRap2EE_PDF = new TH1D("hDDilPtRap2EE_PDF", "hDDilPtRap2EE_PDF", nBinPtRap2, xbinsPtRap2);
  hDDilPtRap3MM_PDF = new TH1D("hDDilPtRap3MM_PDF", "hDDilPtRap3MM_PDF", nBinPtRap3, xbinsPtRap3);
  hDDilPtRap3EE_PDF = new TH1D("hDDilPtRap3EE_PDF", "hDDilPtRap3EE_PDF", nBinPtRap3, xbinsPtRap3);
  hDDilPtRap4MM_PDF = new TH1D("hDDilPtRap4MM_PDF", "hDDilPtRap4MM_PDF", nBinPtRap4, xbinsPtRap4);
  hDDilPtRap4EE_PDF = new TH1D("hDDilPtRap4EE_PDF", "hDDilPtRap4EE_PDF", nBinPtRap4, xbinsPtRap4);
  hDWWMLL_PDF = new TH1D("hDWWMLL_PDF", "hDWWMLL_PDF", nBinWW, xbinsWWMLL);
  hDWWPTL1_PDF = new TH1D("hDWWPTL1_PDF", "hDWWPTL1_PDF", nBinWW, xbinsWWPTL1);
  hDWWPTL2_PDF = new TH1D("hDWWPTL2_PDF", "hDWWPTL2_PDF", nBinWW, xbinsWWPTL2);
  hDWWDPHILL_PDF = new TH1D("hDWWDPHILL_PDF", "hDWWDPHILL_PDF", nBinWW, xbinsWWDPHILL);
  hDWWN0JET_PDF = new TH1D("hDWWN0JET_PDF", "hDWWN0JET_PDF", nBinWWN0JET, xbinsWWN0JET);

  hDDilPtMM_QCD = new TH1D("hDDilPtMM_QCD", "hDDilPtMM_QCD", nBinPt, xbinsPt);
  hDDilPtEE_QCD = new TH1D("hDDilPtEE_QCD", "hDDilPtEE_QCD", nBinPt, xbinsPt);
  hDDilRapMM_QCD = new TH1D("hDDilRapMM_QCD", "hDDilRapMM_QCD", nBinRap, xbinsRap);
  hDDilRapEE_QCD = new TH1D("hDDilRapEE_QCD", "hDDilRapEE_QCD", nBinRap, xbinsRap);
  hDDilPhiStarMM_QCD = new TH1D("hDDilPhiStarMM_QCD", "hDDilPhiStarMM_QCD", nBinPhiStar, xbinsPhiStar);
  hDDilPhiStarEE_QCD = new TH1D("hDDilPhiStarEE_QCD", "hDDilPhiStarEE_QCD", nBinPhiStar, xbinsPhiStar);
  hDDilPtRap0MM_QCD = new TH1D("hDDilPtRap0MM_QCD", "hDDilPtRap0MM_QCD", nBinPtRap0, xbinsPtRap0);
  hDDilPtRap0EE_QCD = new TH1D("hDDilPtRap0EE_QCD", "hDDilPtRap0EE_QCD", nBinPtRap0, xbinsPtRap0);
  hDDilPtRap1MM_QCD = new TH1D("hDDilPtRap1MM_QCD", "hDDilPtRap1MM_QCD", nBinPtRap1, xbinsPtRap1);
  hDDilPtRap1EE_QCD = new TH1D("hDDilPtRap1EE_QCD", "hDDilPtRap1EE_QCD", nBinPtRap1, xbinsPtRap1);
  hDDilPtRap2MM_QCD = new TH1D("hDDilPtRap2MM_QCD", "hDDilPtRap2MM_QCD", nBinPtRap2, xbinsPtRap2);
  hDDilPtRap2EE_QCD = new TH1D("hDDilPtRap2EE_QCD", "hDDilPtRap2EE_QCD", nBinPtRap2, xbinsPtRap2);
  hDDilPtRap3MM_QCD = new TH1D("hDDilPtRap3MM_QCD", "hDDilPtRap3MM_QCD", nBinPtRap3, xbinsPtRap3);
  hDDilPtRap3EE_QCD = new TH1D("hDDilPtRap3EE_QCD", "hDDilPtRap3EE_QCD", nBinPtRap3, xbinsPtRap3);
  hDDilPtRap4MM_QCD = new TH1D("hDDilPtRap4MM_QCD", "hDDilPtRap4MM_QCD", nBinPtRap4, xbinsPtRap4);
  hDDilPtRap4EE_QCD = new TH1D("hDDilPtRap4EE_QCD", "hDDilPtRap4EE_QCD", nBinPtRap4, xbinsPtRap4);
  hDWWMLL_QCD = new TH1D("hDWWMLL_QCD", "hDWWMLL_QCD", nBinWW, xbinsWWMLL);
  hDWWPTL1_QCD = new TH1D("hDWWPTL1_QCD", "hDWWPTL1_QCD", nBinWW, xbinsWWPTL1);
  hDWWPTL2_QCD = new TH1D("hDWWPTL2_QCD", "hDWWPTL2_QCD", nBinWW, xbinsWWPTL2);
  hDWWDPHILL_QCD = new TH1D("hDWWDPHILL_QCD", "hDWWDPHILL_QCD", nBinWW, xbinsWWDPHILL);
  hDWWN0JET_QCD = new TH1D("hDWWN0JET_QCD", "hDWWN0JET_QCD", nBinWWN0JET, xbinsWWN0JET);

  for(int i=0; i<6; i++) hDDilPtMM_QCDPart[i] = new TH1D(Form("hDDilPtMM_QCD_%d",i) ,Form("hDDilPtMM_QCD_%d",i), nBinPt, xbinsPt);
  for(int i=0; i<6; i++) hDDilPtEE_QCDPart[i] = new TH1D(Form("hDDilPtEE_QCD_%d",i) ,Form("hDDilPtEE_QCD_%d",i), nBinPt, xbinsPt);
  for(int i=0; i<6; i++) hDDilRapMM_QCDPart[i] = new TH1D(Form("hDDilRapMM_QCD_%d",i) ,Form("hDDilRapMM_QCD_%d",i), nBinRap, xbinsRap);
  for(int i=0; i<6; i++) hDDilRapEE_QCDPart[i] = new TH1D(Form("hDDilRapEE_QCD_%d",i) ,Form("hDDilRapEE_QCD_%d",i), nBinRap, xbinsRap);
  for(int i=0; i<6; i++) hDDilPhiStarMM_QCDPart[i] = new TH1D(Form("hDDilPhiStarMM_QCD_%d",i) ,Form("hDDilPhiStarMM_QCD_%d",i), nBinPhiStar, xbinsPhiStar);
  for(int i=0; i<6; i++) hDDilPhiStarEE_QCDPart[i] = new TH1D(Form("hDDilPhiStarEE_QCD_%d",i) ,Form("hDDilPhiStarEE_QCD_%d",i), nBinPhiStar, xbinsPhiStar);
  for(int i=0; i<6; i++) hDDilPtRap0MM_QCDPart[i] = new TH1D(Form("hDDilPtRap0MM_QCD_%d",i) ,Form("hDDilPtRap0MM_QCD_%d",i), nBinPtRap0, xbinsPtRap0);
  for(int i=0; i<6; i++) hDDilPtRap0EE_QCDPart[i] = new TH1D(Form("hDDilPtRap0EE_QCD_%d",i) ,Form("hDDilPtRap0EE_QCD_%d",i), nBinPtRap0, xbinsPtRap0);
  for(int i=0; i<6; i++) hDDilPtRap1MM_QCDPart[i] = new TH1D(Form("hDDilPtRap1MM_QCD_%d",i) ,Form("hDDilPtRap1MM_QCD_%d",i), nBinPtRap1, xbinsPtRap1);
  for(int i=0; i<6; i++) hDDilPtRap1EE_QCDPart[i] = new TH1D(Form("hDDilPtRap1EE_QCD_%d",i) ,Form("hDDilPtRap1EE_QCD_%d",i), nBinPtRap1, xbinsPtRap1);
  for(int i=0; i<6; i++) hDDilPtRap2MM_QCDPart[i] = new TH1D(Form("hDDilPtRap2MM_QCD_%d",i) ,Form("hDDilPtRap2MM_QCD_%d",i), nBinPtRap2, xbinsPtRap2);
  for(int i=0; i<6; i++) hDDilPtRap2EE_QCDPart[i] = new TH1D(Form("hDDilPtRap2EE_QCD_%d",i) ,Form("hDDilPtRap2EE_QCD_%d",i), nBinPtRap2, xbinsPtRap2);
  for(int i=0; i<6; i++) hDDilPtRap3MM_QCDPart[i] = new TH1D(Form("hDDilPtRap3MM_QCD_%d",i) ,Form("hDDilPtRap3MM_QCD_%d",i), nBinPtRap3, xbinsPtRap3);
  for(int i=0; i<6; i++) hDDilPtRap3EE_QCDPart[i] = new TH1D(Form("hDDilPtRap3EE_QCD_%d",i) ,Form("hDDilPtRap3EE_QCD_%d",i), nBinPtRap3, xbinsPtRap3);
  for(int i=0; i<6; i++) hDDilPtRap4MM_QCDPart[i] = new TH1D(Form("hDDilPtRap4MM_QCD_%d",i) ,Form("hDDilPtRap4MM_QCD_%d",i), nBinPtRap4, xbinsPtRap4);
  for(int i=0; i<6; i++) hDDilPtRap4EE_QCDPart[i] = new TH1D(Form("hDDilPtRap4EE_QCD_%d",i) ,Form("hDDilPtRap4EE_QCD_%d",i), nBinPtRap4, xbinsPtRap4);
  for(int i=0; i<6; i++) hDWWMLL_QCDPart[i] = new TH1D(Form("hDWWMLL_QCD_%d",i) ,Form("hDWWMLL_QCD_%d",i), nBinWW, xbinsWWMLL);
  for(int i=0; i<6; i++) hDWWPTL1_QCDPart[i] = new TH1D(Form("hDWWPTL1_QCD_%d",i) ,Form("hDWWPTL1_QCD_%d",i), nBinWW, xbinsWWPTL1);
  for(int i=0; i<6; i++) hDWWPTL2_QCDPart[i] = new TH1D(Form("hDWWPTL2_QCD_%d",i) ,Form("hDWWPTL2_QCD_%d",i), nBinWW, xbinsWWPTL2);
  for(int i=0; i<6; i++) hDWWDPHILL_QCDPart[i] = new TH1D(Form("hDWWDPHILL_QCD_%d",i) ,Form("hDWWDPHILL_QCD_%d",i), nBinWW, xbinsWWDPHILL);
  for(int i=0; i<6; i++) hDWWN0JET_QCDPart[i] = new TH1D(Form("hDWWN0JET_QCD_%d",i) ,Form("hDWWN0JET_QCD_%d",i), nBinWWN0JET, xbinsWWN0JET);

  if (weightNames) {
    if (weightNames->GetEntries()!=377 && weightNames->GetEntries()!=22) {
      PError("PandaLeptonicAnalyzer::Init",
          TString::Format("Reweighting failed because only found %u weights!",
                          unsigned(weightNames->GetEntries())));
      return 1;
    }
    TString *id = new TString();
    weightNames->SetBranchAddress("id",&id);
    unsigned nW = weightNames->GetEntriesFast();
    for (unsigned iW=0; iW!=nW; ++iW) {
      weightNames->GetEntry(iW);
      wIDs.push_back(*id);
    }
  } else if (processType==kSignal) {
    PError("PandaLeptonicAnalyzer::Init","This is a signal file, but the weights are missing!");
    return 2;
  }


  // manipulate the output tree
  if (isData) {
    std::vector<TString> droppable = {"mcWeight","scale","pdf.*","gen.*","sf_.*"};
    gt->RemoveBranches(droppable,{"sf_phoPurity"});
  }
  if (flags["lepton"]) {
    std::vector<TString> droppable = {"sf_sjbtag*"};
    gt->RemoveBranches(droppable);
  }
  if (flags["genOnly"]) {
    std::vector<TString> keepable = {"mcWeight","scale","scaleUp",
                                     "scaleDown","pdf*","gen*",
                                     "sf_tt*","sf_qcdTT*"};
    gt->RemoveBranches({".*"},keepable);
  }

  if (DEBUG) PDebug("PandaLeptonicAnalyzer::Init","Finished configuration");

  return 0;
}


panda::GenParticle const *PandaLeptonicAnalyzer::MatchToGen(double eta, double phi, double radius, int pdgid) {
  panda::GenParticle const* found=NULL;
  double r2 = radius*radius;
  pdgid = abs(pdgid);

  unsigned int counter=0;
  for (map<panda::GenParticle const*,float>::iterator iG=genObjects.begin();
      iG!=genObjects.end(); ++iG) {
    if (found!=NULL)
      break;
    if (pdgid!=0 && abs(iG->first->pdgid)!=pdgid)
      continue;
    if (DeltaR2(eta,phi,iG->first->eta(),iG->first->phi())<r2)
      found = iG->first;
  }

  return found;
}


void PandaLeptonicAnalyzer::Terminate() {
  {
    printf("hDDilPtMM: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilPtMM_QCDPart[0]->GetSumOfWeights(),hDDilPtMM_QCDPart[1]->GetSumOfWeights(),hDDilPtMM_QCDPart[2]->GetSumOfWeights(),
  	    hDDilPtMM_QCDPart[3]->GetSumOfWeights(),hDDilPtMM_QCDPart[4]->GetSumOfWeights(),hDDilPtMM_QCDPart[5]->GetSumOfWeights(),hDDilPtMM->GetSumOfWeights());
    for(int nb=1; nb<=nBinPt+1; nb++){
      double systQCDScale = TMath::Abs(hDDilPtMM_QCDPart[0]->GetBinContent(nb)-hDDilPtMM->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilPtMM_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtMM->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilPtMM_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtMM->GetBinContent(nb));
      }
      if(hDDilPtMM->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilPtMM->GetBinContent(nb); else systQCDScale = 1;

      hDDilPtMM_QCD->SetBinContent(nb, hDDilPtMM->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilPtEE: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilPtEE_QCDPart[0]->GetSumOfWeights(),hDDilPtEE_QCDPart[1]->GetSumOfWeights(),hDDilPtEE_QCDPart[2]->GetSumOfWeights(),
  	    hDDilPtEE_QCDPart[3]->GetSumOfWeights(),hDDilPtEE_QCDPart[4]->GetSumOfWeights(),hDDilPtEE_QCDPart[5]->GetSumOfWeights(),hDDilPtEE->GetSumOfWeights());
    for(int nb=1; nb<=nBinPt+1; nb++){
      double systQCDScale = TMath::Abs(hDDilPtEE_QCDPart[0]->GetBinContent(nb)-hDDilPtEE->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilPtEE_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtEE->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilPtEE_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtEE->GetBinContent(nb));
      }
      if(hDDilPtEE->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilPtEE->GetBinContent(nb); else systQCDScale = 1;

      hDDilPtEE_QCD->SetBinContent(nb, hDDilPtEE->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilRapMM: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilRapMM_QCDPart[0]->GetSumOfWeights(),hDDilRapMM_QCDPart[1]->GetSumOfWeights(),hDDilRapMM_QCDPart[2]->GetSumOfWeights(),
  	    hDDilRapMM_QCDPart[3]->GetSumOfWeights(),hDDilRapMM_QCDPart[4]->GetSumOfWeights(),hDDilRapMM_QCDPart[5]->GetSumOfWeights(),hDDilRapMM->GetSumOfWeights());
    for(int nb=1; nb<=nBinRap+1; nb++){
      double systQCDScale = TMath::Abs(hDDilRapMM_QCDPart[0]->GetBinContent(nb)-hDDilRapMM->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilRapMM_QCDPart[nqcd]->GetBinContent(nb)-hDDilRapMM->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilRapMM_QCDPart[nqcd]->GetBinContent(nb)-hDDilRapMM->GetBinContent(nb));
      }
      if(hDDilRapMM->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilRapMM->GetBinContent(nb); else systQCDScale = 1;

      hDDilRapMM_QCD->SetBinContent(nb, hDDilRapMM->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilRapEE: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilRapEE_QCDPart[0]->GetSumOfWeights(),hDDilRapEE_QCDPart[1]->GetSumOfWeights(),hDDilRapEE_QCDPart[2]->GetSumOfWeights(),
  	    hDDilRapEE_QCDPart[3]->GetSumOfWeights(),hDDilRapEE_QCDPart[4]->GetSumOfWeights(),hDDilRapEE_QCDPart[5]->GetSumOfWeights(),hDDilRapEE->GetSumOfWeights());
    for(int nb=1; nb<=nBinRap+1; nb++){
      double systQCDScale = TMath::Abs(hDDilRapEE_QCDPart[0]->GetBinContent(nb)-hDDilRapEE->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilRapEE_QCDPart[nqcd]->GetBinContent(nb)-hDDilRapEE->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilRapEE_QCDPart[nqcd]->GetBinContent(nb)-hDDilRapEE->GetBinContent(nb));
      }
      if(hDDilRapEE->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilRapEE->GetBinContent(nb); else systQCDScale = 1;

      hDDilRapEE_QCD->SetBinContent(nb, hDDilRapEE->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilPhiStarMM: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilPhiStarMM_QCDPart[0]->GetSumOfWeights(),hDDilPhiStarMM_QCDPart[1]->GetSumOfWeights(),hDDilPhiStarMM_QCDPart[2]->GetSumOfWeights(),
  	    hDDilPhiStarMM_QCDPart[3]->GetSumOfWeights(),hDDilPhiStarMM_QCDPart[4]->GetSumOfWeights(),hDDilPhiStarMM_QCDPart[5]->GetSumOfWeights(),hDDilPhiStarMM->GetSumOfWeights());
    for(int nb=1; nb<=nBinPhiStar+1; nb++){
      double systQCDScale = TMath::Abs(hDDilPhiStarMM_QCDPart[0]->GetBinContent(nb)-hDDilPhiStarMM->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilPhiStarMM_QCDPart[nqcd]->GetBinContent(nb)-hDDilPhiStarMM->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilPhiStarMM_QCDPart[nqcd]->GetBinContent(nb)-hDDilPhiStarMM->GetBinContent(nb));
      }
      if(hDDilPhiStarMM->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilPhiStarMM->GetBinContent(nb); else systQCDScale = 1;

      hDDilPhiStarMM_QCD->SetBinContent(nb, hDDilPhiStarMM->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilPhiStarEE: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilPhiStarEE_QCDPart[0]->GetSumOfWeights(),hDDilPhiStarEE_QCDPart[1]->GetSumOfWeights(),hDDilPhiStarEE_QCDPart[2]->GetSumOfWeights(),
  	    hDDilPhiStarEE_QCDPart[3]->GetSumOfWeights(),hDDilPhiStarEE_QCDPart[4]->GetSumOfWeights(),hDDilPhiStarEE_QCDPart[5]->GetSumOfWeights(),hDDilPhiStarEE->GetSumOfWeights());
    for(int nb=1; nb<=nBinPhiStar+1; nb++){
      double systQCDScale = TMath::Abs(hDDilPhiStarEE_QCDPart[0]->GetBinContent(nb)-hDDilPhiStarEE->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilPhiStarEE_QCDPart[nqcd]->GetBinContent(nb)-hDDilPhiStarEE->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilPhiStarEE_QCDPart[nqcd]->GetBinContent(nb)-hDDilPhiStarEE->GetBinContent(nb));
      }
      if(hDDilPhiStarEE->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilPhiStarEE->GetBinContent(nb); else systQCDScale = 1;

      hDDilPhiStarEE_QCD->SetBinContent(nb, hDDilPhiStarEE->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilPtRap0MM: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilPtRap0MM_QCDPart[0]->GetSumOfWeights(),hDDilPtRap0MM_QCDPart[1]->GetSumOfWeights(),hDDilPtRap0MM_QCDPart[2]->GetSumOfWeights(),
  	    hDDilPtRap0MM_QCDPart[3]->GetSumOfWeights(),hDDilPtRap0MM_QCDPart[4]->GetSumOfWeights(),hDDilPtRap0MM_QCDPart[5]->GetSumOfWeights(),hDDilPtRap0MM->GetSumOfWeights());
    for(int nb=1; nb<=nBinPtRap0+1; nb++){
      double systQCDScale = TMath::Abs(hDDilPtRap0MM_QCDPart[0]->GetBinContent(nb)-hDDilPtRap0MM->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilPtRap0MM_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap0MM->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilPtRap0MM_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap0MM->GetBinContent(nb));
      }
      if(hDDilPtRap0MM->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilPtRap0MM->GetBinContent(nb); else systQCDScale = 1;

      hDDilPtRap0MM_QCD->SetBinContent(nb, hDDilPtRap0MM->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilPtRap0EE: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilPtRap0EE_QCDPart[0]->GetSumOfWeights(),hDDilPtRap0EE_QCDPart[1]->GetSumOfWeights(),hDDilPtRap0EE_QCDPart[2]->GetSumOfWeights(),
  	    hDDilPtRap0EE_QCDPart[3]->GetSumOfWeights(),hDDilPtRap0EE_QCDPart[4]->GetSumOfWeights(),hDDilPtRap0EE_QCDPart[5]->GetSumOfWeights(),hDDilPtRap0EE->GetSumOfWeights());
    for(int nb=1; nb<=nBinPtRap0+1; nb++){
      double systQCDScale = TMath::Abs(hDDilPtRap0EE_QCDPart[0]->GetBinContent(nb)-hDDilPtRap0EE->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilPtRap0EE_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap0EE->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilPtRap0EE_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap0EE->GetBinContent(nb));
      }
      if(hDDilPtRap0EE->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilPtRap0EE->GetBinContent(nb); else systQCDScale = 1;

      hDDilPtRap0EE_QCD->SetBinContent(nb, hDDilPtRap0EE->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilPtRap1MM: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilPtRap1MM_QCDPart[0]->GetSumOfWeights(),hDDilPtRap1MM_QCDPart[1]->GetSumOfWeights(),hDDilPtRap1MM_QCDPart[2]->GetSumOfWeights(),
  	    hDDilPtRap1MM_QCDPart[3]->GetSumOfWeights(),hDDilPtRap1MM_QCDPart[4]->GetSumOfWeights(),hDDilPtRap1MM_QCDPart[5]->GetSumOfWeights(),hDDilPtRap1MM->GetSumOfWeights());
    for(int nb=1; nb<=nBinPtRap1+1; nb++){
      double systQCDScale = TMath::Abs(hDDilPtRap1MM_QCDPart[0]->GetBinContent(nb)-hDDilPtRap1MM->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilPtRap1MM_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap1MM->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilPtRap1MM_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap1MM->GetBinContent(nb));
      }
      if(hDDilPtRap1MM->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilPtRap1MM->GetBinContent(nb); else systQCDScale = 1;

      hDDilPtRap1MM_QCD->SetBinContent(nb, hDDilPtRap1MM->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilPtRap1EE: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilPtRap1EE_QCDPart[0]->GetSumOfWeights(),hDDilPtRap1EE_QCDPart[1]->GetSumOfWeights(),hDDilPtRap1EE_QCDPart[2]->GetSumOfWeights(),
  	    hDDilPtRap1EE_QCDPart[3]->GetSumOfWeights(),hDDilPtRap1EE_QCDPart[4]->GetSumOfWeights(),hDDilPtRap1EE_QCDPart[5]->GetSumOfWeights(),hDDilPtRap1EE->GetSumOfWeights());
    for(int nb=1; nb<=nBinPtRap1+1; nb++){
      double systQCDScale = TMath::Abs(hDDilPtRap1EE_QCDPart[0]->GetBinContent(nb)-hDDilPtRap1EE->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilPtRap1EE_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap1EE->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilPtRap1EE_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap1EE->GetBinContent(nb));
      }
      if(hDDilPtRap1EE->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilPtRap1EE->GetBinContent(nb); else systQCDScale = 1;

      hDDilPtRap1EE_QCD->SetBinContent(nb, hDDilPtRap1EE->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilPtRap2MM: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilPtRap2MM_QCDPart[0]->GetSumOfWeights(),hDDilPtRap2MM_QCDPart[1]->GetSumOfWeights(),hDDilPtRap2MM_QCDPart[2]->GetSumOfWeights(),
  	    hDDilPtRap2MM_QCDPart[3]->GetSumOfWeights(),hDDilPtRap2MM_QCDPart[4]->GetSumOfWeights(),hDDilPtRap2MM_QCDPart[5]->GetSumOfWeights(),hDDilPtRap2MM->GetSumOfWeights());
    for(int nb=1; nb<=nBinPtRap2+1; nb++){
      double systQCDScale = TMath::Abs(hDDilPtRap2MM_QCDPart[0]->GetBinContent(nb)-hDDilPtRap2MM->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilPtRap2MM_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap2MM->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilPtRap2MM_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap2MM->GetBinContent(nb));
      }
      if(hDDilPtRap2MM->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilPtRap2MM->GetBinContent(nb); else systQCDScale = 1;

      hDDilPtRap2MM_QCD->SetBinContent(nb, hDDilPtRap2MM->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilPtRap2EE: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilPtRap2EE_QCDPart[0]->GetSumOfWeights(),hDDilPtRap2EE_QCDPart[1]->GetSumOfWeights(),hDDilPtRap2EE_QCDPart[2]->GetSumOfWeights(),
  	    hDDilPtRap2EE_QCDPart[3]->GetSumOfWeights(),hDDilPtRap2EE_QCDPart[4]->GetSumOfWeights(),hDDilPtRap2EE_QCDPart[5]->GetSumOfWeights(),hDDilPtRap2EE->GetSumOfWeights());
    for(int nb=1; nb<=nBinPtRap2+1; nb++){
      double systQCDScale = TMath::Abs(hDDilPtRap2EE_QCDPart[0]->GetBinContent(nb)-hDDilPtRap2EE->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilPtRap2EE_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap2EE->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilPtRap2EE_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap2EE->GetBinContent(nb));
      }
      if(hDDilPtRap2EE->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilPtRap2EE->GetBinContent(nb); else systQCDScale = 1;

      hDDilPtRap2EE_QCD->SetBinContent(nb, hDDilPtRap2EE->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilPtRap3MM: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilPtRap3MM_QCDPart[0]->GetSumOfWeights(),hDDilPtRap3MM_QCDPart[1]->GetSumOfWeights(),hDDilPtRap3MM_QCDPart[2]->GetSumOfWeights(),
  	    hDDilPtRap3MM_QCDPart[3]->GetSumOfWeights(),hDDilPtRap3MM_QCDPart[4]->GetSumOfWeights(),hDDilPtRap3MM_QCDPart[5]->GetSumOfWeights(),hDDilPtRap3MM->GetSumOfWeights());
    for(int nb=1; nb<=nBinPtRap3+1; nb++){
      double systQCDScale = TMath::Abs(hDDilPtRap3MM_QCDPart[0]->GetBinContent(nb)-hDDilPtRap3MM->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilPtRap3MM_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap3MM->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilPtRap3MM_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap3MM->GetBinContent(nb));
      }
      if(hDDilPtRap3MM->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilPtRap3MM->GetBinContent(nb); else systQCDScale = 1;

      hDDilPtRap3MM_QCD->SetBinContent(nb, hDDilPtRap3MM->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilPtRap3EE: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilPtRap3EE_QCDPart[0]->GetSumOfWeights(),hDDilPtRap3EE_QCDPart[1]->GetSumOfWeights(),hDDilPtRap3EE_QCDPart[2]->GetSumOfWeights(),
  	    hDDilPtRap3EE_QCDPart[3]->GetSumOfWeights(),hDDilPtRap3EE_QCDPart[4]->GetSumOfWeights(),hDDilPtRap3EE_QCDPart[5]->GetSumOfWeights(),hDDilPtRap3EE->GetSumOfWeights());
    for(int nb=1; nb<=nBinPtRap3+1; nb++){
      double systQCDScale = TMath::Abs(hDDilPtRap3EE_QCDPart[0]->GetBinContent(nb)-hDDilPtRap3EE->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilPtRap3EE_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap3EE->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilPtRap3EE_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap3EE->GetBinContent(nb));
      }
      if(hDDilPtRap3EE->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilPtRap3EE->GetBinContent(nb); else systQCDScale = 1;

      hDDilPtRap3EE_QCD->SetBinContent(nb, hDDilPtRap3EE->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilPtRap4MM: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilPtRap4MM_QCDPart[0]->GetSumOfWeights(),hDDilPtRap4MM_QCDPart[1]->GetSumOfWeights(),hDDilPtRap4MM_QCDPart[2]->GetSumOfWeights(),
  	    hDDilPtRap4MM_QCDPart[3]->GetSumOfWeights(),hDDilPtRap4MM_QCDPart[4]->GetSumOfWeights(),hDDilPtRap4MM_QCDPart[5]->GetSumOfWeights(),hDDilPtRap4MM->GetSumOfWeights());
    for(int nb=1; nb<=nBinPtRap4+1; nb++){
      double systQCDScale = TMath::Abs(hDDilPtRap4MM_QCDPart[0]->GetBinContent(nb)-hDDilPtRap4MM->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilPtRap4MM_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap4MM->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilPtRap4MM_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap4MM->GetBinContent(nb));
      }
      if(hDDilPtRap4MM->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilPtRap4MM->GetBinContent(nb); else systQCDScale = 1;

      hDDilPtRap4MM_QCD->SetBinContent(nb, hDDilPtRap4MM->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDDilPtRap4EE: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDDilPtRap4EE_QCDPart[0]->GetSumOfWeights(),hDDilPtRap4EE_QCDPart[1]->GetSumOfWeights(),hDDilPtRap4EE_QCDPart[2]->GetSumOfWeights(),
  	    hDDilPtRap4EE_QCDPart[3]->GetSumOfWeights(),hDDilPtRap4EE_QCDPart[4]->GetSumOfWeights(),hDDilPtRap4EE_QCDPart[5]->GetSumOfWeights(),hDDilPtRap4EE->GetSumOfWeights());
    for(int nb=1; nb<=nBinPtRap4+1; nb++){
      double systQCDScale = TMath::Abs(hDDilPtRap4EE_QCDPart[0]->GetBinContent(nb)-hDDilPtRap4EE->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDDilPtRap4EE_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap4EE->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDDilPtRap4EE_QCDPart[nqcd]->GetBinContent(nb)-hDDilPtRap4EE->GetBinContent(nb));
      }
      if(hDDilPtRap4EE->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDDilPtRap4EE->GetBinContent(nb); else systQCDScale = 1;

      hDDilPtRap4EE_QCD->SetBinContent(nb, hDDilPtRap4EE->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDWWSSMLL: %f\n",hDWWSSMLL->GetSumOfWeights());
  }
  {
    printf("hDWWMLL: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDWWMLL_QCDPart[0]->GetSumOfWeights(),hDWWMLL_QCDPart[1]->GetSumOfWeights(),hDWWMLL_QCDPart[2]->GetSumOfWeights(),
  	    hDWWMLL_QCDPart[3]->GetSumOfWeights(),hDWWMLL_QCDPart[4]->GetSumOfWeights(),hDWWMLL_QCDPart[5]->GetSumOfWeights(),hDWWMLL->GetSumOfWeights());
    for(int nb=1; nb<=nBinPt+1; nb++){
      double systQCDScale = TMath::Abs(hDWWMLL_QCDPart[0]->GetBinContent(nb)-hDWWMLL->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDWWMLL_QCDPart[nqcd]->GetBinContent(nb)-hDWWMLL->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDWWMLL_QCDPart[nqcd]->GetBinContent(nb)-hDWWMLL->GetBinContent(nb));
      }
      if(hDWWMLL->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDWWMLL->GetBinContent(nb); else systQCDScale = 1;

      hDWWMLL_QCD->SetBinContent(nb, hDWWMLL->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDWWPTL1: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDWWPTL1_QCDPart[0]->GetSumOfWeights(),hDWWPTL1_QCDPart[1]->GetSumOfWeights(),hDWWPTL1_QCDPart[2]->GetSumOfWeights(),
  	    hDWWPTL1_QCDPart[3]->GetSumOfWeights(),hDWWPTL1_QCDPart[4]->GetSumOfWeights(),hDWWPTL1_QCDPart[5]->GetSumOfWeights(),hDWWPTL1->GetSumOfWeights());
    for(int nb=1; nb<=nBinPt+1; nb++){
      double systQCDScale = TMath::Abs(hDWWPTL1_QCDPart[0]->GetBinContent(nb)-hDWWPTL1->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDWWPTL1_QCDPart[nqcd]->GetBinContent(nb)-hDWWPTL1->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDWWPTL1_QCDPart[nqcd]->GetBinContent(nb)-hDWWPTL1->GetBinContent(nb));
      }
      if(hDWWPTL1->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDWWPTL1->GetBinContent(nb); else systQCDScale = 1;

      hDWWPTL1_QCD->SetBinContent(nb, hDWWPTL1->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDWWPTL2: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDWWPTL2_QCDPart[0]->GetSumOfWeights(),hDWWPTL2_QCDPart[1]->GetSumOfWeights(),hDWWPTL2_QCDPart[2]->GetSumOfWeights(),
  	    hDWWPTL2_QCDPart[3]->GetSumOfWeights(),hDWWPTL2_QCDPart[4]->GetSumOfWeights(),hDWWPTL2_QCDPart[5]->GetSumOfWeights(),hDWWPTL2->GetSumOfWeights());
    for(int nb=1; nb<=nBinPt+1; nb++){
      double systQCDScale = TMath::Abs(hDWWPTL2_QCDPart[0]->GetBinContent(nb)-hDWWPTL2->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDWWPTL2_QCDPart[nqcd]->GetBinContent(nb)-hDWWPTL2->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDWWPTL2_QCDPart[nqcd]->GetBinContent(nb)-hDWWPTL2->GetBinContent(nb));
      }
      if(hDWWPTL2->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDWWPTL2->GetBinContent(nb); else systQCDScale = 1;

      hDWWPTL2_QCD->SetBinContent(nb, hDWWPTL2->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDWWDPHILL: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDWWDPHILL_QCDPart[0]->GetSumOfWeights(),hDWWDPHILL_QCDPart[1]->GetSumOfWeights(),hDWWDPHILL_QCDPart[2]->GetSumOfWeights(),
  	    hDWWDPHILL_QCDPart[3]->GetSumOfWeights(),hDWWDPHILL_QCDPart[4]->GetSumOfWeights(),hDWWDPHILL_QCDPart[5]->GetSumOfWeights(),hDWWDPHILL->GetSumOfWeights());
    for(int nb=1; nb<=nBinPt+1; nb++){
      double systQCDScale = TMath::Abs(hDWWDPHILL_QCDPart[0]->GetBinContent(nb)-hDWWDPHILL->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDWWDPHILL_QCDPart[nqcd]->GetBinContent(nb)-hDWWDPHILL->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDWWDPHILL_QCDPart[nqcd]->GetBinContent(nb)-hDWWDPHILL->GetBinContent(nb));
      }
      if(hDWWDPHILL->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDWWDPHILL->GetBinContent(nb); else systQCDScale = 1;

      hDWWDPHILL_QCD->SetBinContent(nb, hDWWDPHILL->GetBinContent(nb)*systQCDScale);
    }
  }
  {
    printf("hDWWN0JET: (%f/%f/%f/%f/%f/%f->%f)\n",
  	    hDWWN0JET_QCDPart[0]->GetSumOfWeights(),hDWWN0JET_QCDPart[1]->GetSumOfWeights(),hDWWN0JET_QCDPart[2]->GetSumOfWeights(),
  	    hDWWN0JET_QCDPart[3]->GetSumOfWeights(),hDWWN0JET_QCDPart[4]->GetSumOfWeights(),hDWWN0JET_QCDPart[5]->GetSumOfWeights(),hDWWN0JET->GetSumOfWeights());
    for(int nb=1; nb<=nBinPt+1; nb++){
      double systQCDScale = TMath::Abs(hDWWN0JET_QCDPart[0]->GetBinContent(nb)-hDWWN0JET->GetBinContent(nb));

      for(int nqcd=1; nqcd<6; nqcd++) {
        if(TMath::Abs(hDWWN0JET_QCDPart[nqcd]->GetBinContent(nb)-hDWWN0JET->GetBinContent(nb)) > systQCDScale) systQCDScale = TMath::Abs(hDWWN0JET_QCDPart[nqcd]->GetBinContent(nb)-hDWWN0JET->GetBinContent(nb));
      }
      if(hDWWN0JET->GetBinContent(nb) > 0) systQCDScale = 1.0+systQCDScale/hDWWN0JET->GetBinContent(nb); else systQCDScale = 1;

      hDWWN0JET_QCD->SetBinContent(nb, hDWWN0JET->GetBinContent(nb)*systQCDScale);
    }
  }



  fOut->WriteTObject(tOut);
  fOut->WriteTObject(hDDilPtMM);     fOut->WriteTObject(hDDilPtMM_PDF);     fOut->WriteTObject(hDDilPtMM_QCD);    
  fOut->WriteTObject(hDDilPtEE);     fOut->WriteTObject(hDDilPtEE_PDF);     fOut->WriteTObject(hDDilPtEE_QCD);    
  fOut->WriteTObject(hDDilRapMM);    fOut->WriteTObject(hDDilRapMM_PDF);    fOut->WriteTObject(hDDilRapMM_QCD);   
  fOut->WriteTObject(hDDilRapEE);    fOut->WriteTObject(hDDilRapEE_PDF);    fOut->WriteTObject(hDDilRapEE_QCD);   
  fOut->WriteTObject(hDDilPhiStarMM);fOut->WriteTObject(hDDilPhiStarMM_PDF);fOut->WriteTObject(hDDilPhiStarMM_QCD);    
  fOut->WriteTObject(hDDilPhiStarEE);fOut->WriteTObject(hDDilPhiStarEE_PDF);fOut->WriteTObject(hDDilPhiStarEE_QCD);    
  fOut->WriteTObject(hDDilPtRap0MM); fOut->WriteTObject(hDDilPtRap0MM_PDF); fOut->WriteTObject(hDDilPtRap0MM_QCD);
  fOut->WriteTObject(hDDilPtRap0EE); fOut->WriteTObject(hDDilPtRap0EE_PDF); fOut->WriteTObject(hDDilPtRap0EE_QCD);
  fOut->WriteTObject(hDDilPtRap1MM); fOut->WriteTObject(hDDilPtRap1MM_PDF); fOut->WriteTObject(hDDilPtRap1MM_QCD);
  fOut->WriteTObject(hDDilPtRap1EE); fOut->WriteTObject(hDDilPtRap1EE_PDF); fOut->WriteTObject(hDDilPtRap1EE_QCD);
  fOut->WriteTObject(hDDilPtRap2MM); fOut->WriteTObject(hDDilPtRap2MM_PDF); fOut->WriteTObject(hDDilPtRap2MM_QCD);
  fOut->WriteTObject(hDDilPtRap2EE); fOut->WriteTObject(hDDilPtRap2EE_PDF); fOut->WriteTObject(hDDilPtRap2EE_QCD);
  fOut->WriteTObject(hDDilPtRap3MM); fOut->WriteTObject(hDDilPtRap3MM_PDF); fOut->WriteTObject(hDDilPtRap3MM_QCD);
  fOut->WriteTObject(hDDilPtRap3EE); fOut->WriteTObject(hDDilPtRap3EE_PDF); fOut->WriteTObject(hDDilPtRap3EE_QCD);
  fOut->WriteTObject(hDDilPtRap4MM); fOut->WriteTObject(hDDilPtRap4MM_PDF); fOut->WriteTObject(hDDilPtRap4MM_QCD);
  fOut->WriteTObject(hDDilPtRap4EE); fOut->WriteTObject(hDDilPtRap4EE_PDF); fOut->WriteTObject(hDDilPtRap4EE_QCD);
  fOut->WriteTObject(hDWWSSMLL);  
  fOut->WriteTObject(hDWWMLL);       fOut->WriteTObject(hDWWMLL_PDF);       fOut->WriteTObject(hDWWMLL_QCD);    
  fOut->WriteTObject(hDWWPTL1);      fOut->WriteTObject(hDWWPTL1_PDF);      fOut->WriteTObject(hDWWPTL1_QCD);    
  fOut->WriteTObject(hDWWPTL2);      fOut->WriteTObject(hDWWPTL2_PDF);      fOut->WriteTObject(hDWWPTL2_QCD);    
  fOut->WriteTObject(hDWWDPHILL);    fOut->WriteTObject(hDWWDPHILL_PDF);    fOut->WriteTObject(hDWWDPHILL_QCD);    
  fOut->WriteTObject(hDWWN0JET);     fOut->WriteTObject(hDWWN0JET_PDF);     fOut->WriteTObject(hDWWN0JET_QCD);    
  fOut->Close();

  //for (auto *f : fCorrs)
  //  if (f)
  //    f->Close();
  for (auto *h : h1Corrs)
    delete h;
  for (auto *h : h2Corrs)
    delete h;

  delete btagCalib;
  for (auto *reader : btagReaders )
    delete reader;

  for (auto& iter : ak8UncReader)
    delete iter.second;

  delete ak8JERReader;

  for (auto& iter : ak4UncReader)
    delete iter.second;

  for (auto& iter : ak4ScaleReader)
    delete iter.second;

  delete ak4JERReader;
  
  delete hDTotalMCWeight;
  delete hDDilPtMM;     delete hDDilPtMM_PDF;	  delete hDDilPtMM_QCD;	    for(int i=0; i<6; i++) delete hDDilPtMM_QCDPart[i];
  delete hDDilPtEE;     delete hDDilPtEE_PDF;	  delete hDDilPtEE_QCD;	    for(int i=0; i<6; i++) delete hDDilPtEE_QCDPart[i];
  delete hDDilRapMM;    delete hDDilRapMM_PDF;    delete hDDilRapMM_QCD;    for(int i=0; i<6; i++) delete hDDilRapMM_QCDPart[i];   
  delete hDDilRapEE;    delete hDDilRapEE_PDF;    delete hDDilRapEE_QCD;    for(int i=0; i<6; i++) delete hDDilRapEE_QCDPart[i];   
  delete hDDilPhiStarMM;delete hDDilPhiStarMM_PDF;delete hDDilPhiStarMM_QCD;for(int i=0; i<6; i++) delete hDDilPhiStarMM_QCDPart[i];	
  delete hDDilPhiStarEE;delete hDDilPhiStarEE_PDF;delete hDDilPhiStarEE_QCD;for(int i=0; i<6; i++) delete hDDilPhiStarEE_QCDPart[i];	
  delete hDDilPtRap0MM; delete hDDilPtRap0MM_PDF; delete hDDilPtRap0MM_QCD; for(int i=0; i<6; i++) delete hDDilPtRap0MM_QCDPart[i];
  delete hDDilPtRap0EE; delete hDDilPtRap0EE_PDF; delete hDDilPtRap0EE_QCD; for(int i=0; i<6; i++) delete hDDilPtRap0EE_QCDPart[i];
  delete hDDilPtRap1MM; delete hDDilPtRap1MM_PDF; delete hDDilPtRap1MM_QCD; for(int i=0; i<6; i++) delete hDDilPtRap1MM_QCDPart[i];
  delete hDDilPtRap1EE; delete hDDilPtRap1EE_PDF; delete hDDilPtRap1EE_QCD; for(int i=0; i<6; i++) delete hDDilPtRap1EE_QCDPart[i];
  delete hDDilPtRap2MM; delete hDDilPtRap2MM_PDF; delete hDDilPtRap2MM_QCD; for(int i=0; i<6; i++) delete hDDilPtRap2MM_QCDPart[i];
  delete hDDilPtRap2EE; delete hDDilPtRap2EE_PDF; delete hDDilPtRap2EE_QCD; for(int i=0; i<6; i++) delete hDDilPtRap2EE_QCDPart[i];
  delete hDDilPtRap3MM; delete hDDilPtRap3MM_PDF; delete hDDilPtRap3MM_QCD; for(int i=0; i<6; i++) delete hDDilPtRap3MM_QCDPart[i];
  delete hDDilPtRap3EE; delete hDDilPtRap3EE_PDF; delete hDDilPtRap3EE_QCD; for(int i=0; i<6; i++) delete hDDilPtRap3EE_QCDPart[i];
  delete hDDilPtRap4MM; delete hDDilPtRap4MM_PDF; delete hDDilPtRap4MM_QCD; for(int i=0; i<6; i++) delete hDDilPtRap4MM_QCDPart[i];
  delete hDDilPtRap4EE; delete hDDilPtRap4EE_PDF; delete hDDilPtRap4EE_QCD; for(int i=0; i<6; i++) delete hDDilPtRap4EE_QCDPart[i];

  if (DEBUG) PDebug("PandaLeptonicAnalyzer::Terminate","Finished with output");
}

void PandaLeptonicAnalyzer::OpenCorrection(CorrectionType ct, TString fpath, TString hname, int dim) {
  fCorrs[ct] = TFile::Open(fpath);
  if (dim==1) 
    h1Corrs[ct] = new THCorr1((TH1D*)fCorrs[ct]->Get(hname));
  else
    h2Corrs[ct] = new THCorr2((TH2D*)fCorrs[ct]->Get(hname));
}

double PandaLeptonicAnalyzer::GetCorr(CorrectionType ct, double x, double y) {
  if (h1Corrs[ct]!=0) {
    return h1Corrs[ct]->Eval(x); 
  } else if (h2Corrs[ct]!=0) {
    return h2Corrs[ct]->Eval(x,y);
  } else {
    PError("PandaLeptonicAnalyzer::GetCorr",
       TString::Format("No correction is defined for CorrectionType=%u",ct));
    return 1;
  }
}

double PandaLeptonicAnalyzer::GetError(CorrectionType ct, double x, double y) {
  if (h1Corrs[ct]!=0) {
    return h1Corrs[ct]->Error(x); 
  } else if (h2Corrs[ct]!=0) {
    return h2Corrs[ct]->Error(x,y);
  } else {
    PError("PandaLeptonicAnalyzer::GetCorr",
       TString::Format("No correction is defined for CorrectionType=%u",ct));
    return 1;
  }
}

void PandaLeptonicAnalyzer::SetDataDir(const char *s2) {
  TString dirPath1 = TString(gSystem->Getenv("CMSSW_BASE")) + "/src/";
  TString dirPath2(s2);

  if (DEBUG) PDebug("PandaLeptonicAnalyzer::SetDataDir","Starting loading of data");

  // pileup
  OpenCorrection(cPU    ,dirPath1+"MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root","puWeights",1);
  OpenCorrection(cPUUp  ,dirPath1+"MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root","puWeightsUp",1);
  OpenCorrection(cPUDown,dirPath1+"MitAnalysisRunII/data/80x/puWeights_80x_37ifb.root","puWeightsDown",1);

  OpenCorrection(ZHEwkCorr    ,dirPath1+"MitAnalysisRunII/data/80x/Zll_nloEWK_weight_unnormalized.root","SignalWeight_nloEWK_rebin",1);
  OpenCorrection(ZHEwkCorrUp  ,dirPath1+"MitAnalysisRunII/data/80x/Zll_nloEWK_weight_unnormalized.root","SignalWeight_nloEWK_up_rebin",1);
  OpenCorrection(ZHEwkCorrDown,dirPath1+"MitAnalysisRunII/data/80x/Zll_nloEWK_weight_unnormalized.root","SignalWeight_nloEWK_down_rebin",1);

  OpenCorrection(cLooseMuonId  ,dirPath1+"MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root","scalefactors_MuonLooseId_Muon",2);

  OpenCorrection(cMediumMuonId ,dirPath1+"MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root","scalefactors_MuonMediumId_Muon",2);

  OpenCorrection(cTightMuonId  ,dirPath1+"MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root","scalefactors_TightId_Muon",2);
  OpenCorrection(cLooseMuonIso ,dirPath1+"MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root","scalefactors_Iso_MuonLooseId",2);
  OpenCorrection(cMediumMuonIso,dirPath1+"MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root","scalefactors_Iso_MuonMediumId",2);
  OpenCorrection(cTightMuonIso ,dirPath1+"MitAnalysisRunII/data/80x/muon_scalefactors_37ifb.root","scalefactors_Iso_MuonTightId",2);
  OpenCorrection(cTrackingMuon ,dirPath1+"MitAnalysisRunII/data/80x/Tracking_EfficienciesAndSF_BCDEFGH.root","ratio_eff_eta3_dr030e030_corr",1);

  OpenCorrection(cLooseElectronId ,dirPath1+"MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root","scalefactors_Loose_Electron",2);

  OpenCorrection(cMediumElectronId,dirPath1+"MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root","scalefactors_Medium_Electron",2);

  OpenCorrection(cTightElectronId ,dirPath1+"MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root","scalefactors_Tight_Electron",2);
  OpenCorrection(cTrackingElectron,dirPath1+"MitAnalysisRunII/data/80x/scalefactors_80x_egpog_37ifb.root","scalefactors_Reco_Electron",2);

  // btag SFs
  btagCalib = new BTagCalibration("csvv2",(dirPath1+"MitAnalysisRunII/data/80x/CSVv2_Moriond17_B_H.csv").Data());
  btagReaders[bJetL] = new BTagCalibrationReader(BTagEntry::OP_LOOSE,"central",{"up","down"});
  btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_B,"comb");
  btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_C,"comb");
  btagReaders[bJetL]->load(*btagCalib,BTagEntry::FLAV_UDSG,"incl");

  btagReaders[bJetM] = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"});
  btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_B,"comb");
  btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_C,"comb");
  btagReaders[bJetM]->load(*btagCalib,BTagEntry::FLAV_UDSG,"incl");

  if (DEBUG) PDebug("PandaLeptonicAnalyzer::SetDataDir","Loaded btag SFs");

  TString jecV = "V4", jecReco = "23Sep2016"; 
  TString jecVFull = jecReco+jecV;
  ak8UncReader["MC"] = new JetCorrectionUncertainty(
     (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_Uncertainty_AK8PFPuppi.txt").Data()
    );
  std::vector<TString> eraGroups = {"BCD","EF","G","H"};
  for (auto e : eraGroups) {
    ak8UncReader["data"+e] = new JetCorrectionUncertainty(
       (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_Uncertainty_AK8PFPuppi.txt").Data()
      );
  }

  ak8JERReader = new JERReader(dirPath2+"/jec/25nsV10/Spring16_25nsV10_MC_SF_AK8PFPuppi.txt",
                               dirPath2+"/jec/25nsV10/Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt");


  ak4UncReader["MC"] = new JetCorrectionUncertainty(
     (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_Uncertainty_AK4PFPuppi.txt").Data()
    );
  for (auto e : eraGroups) {
    ak4UncReader["data"+e] = new JetCorrectionUncertainty(
       (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_Uncertainty_AK4PFPuppi.txt").Data()
      );
  }

  ak4JERReader = new JERReader(dirPath2+"/jec/25nsV10/Spring16_25nsV10_MC_SF_AK4PFPuppi.txt",
                               dirPath2+"/jec/25nsV10/Spring16_25nsV10_MC_PtResolution_AK4PFPuppi.txt");

  std::vector<JetCorrectorParameters> params = {
    JetCorrectorParameters(
      (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L1FastJet_AK4PFPuppi.txt").Data()),
    JetCorrectorParameters(
      (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L2Relative_AK4PFPuppi.txt").Data()),
    JetCorrectorParameters(
      (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L3Absolute_AK4PFPuppi.txt").Data()),
    JetCorrectorParameters(
      (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L2L3Residual_AK4PFPuppi.txt").Data())
  };
  ak4ScaleReader["MC"] = new FactorizedJetCorrector(params);
  if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::SetDataDir","Loaded JES for AK4 MC");
  for (auto e : eraGroups) {
    params = {
      JetCorrectorParameters(
        (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L1FastJet_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L2Relative_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L3Absolute_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath2+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L2L3Residual_AK4PFPuppi.txt").Data())
    };
    ak4ScaleReader["data"+e] = new FactorizedJetCorrector(params);
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::SetDataDir","Loaded JES for AK4 "+e);
  }

  if (DEBUG) PDebug("PandaLeptonicAnalyzer::SetDataDir","Loaded JES/R");

}

void PandaLeptonicAnalyzer::AddGoodLumiRange(int run, int l0, int l1) {
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) { // don't know about this run yet
    std::vector<LumiRange> newLumiList;
    newLumiList.emplace_back(l0,l1);
    goodLumis[run] = newLumiList;
  } else {
    run_->second.emplace_back(l0,l1);
  }
}

bool PandaLeptonicAnalyzer::PassGoodLumis(int run, int lumi) {
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) {
    // matched no run
    if (DEBUG) 
      PDebug("PandaLeptonicAnalyzer::PassGoodLumis",TString::Format("Failing run=%i",run));
    return false;
  }

  // found the run, now look for a lumi range
  for (auto &range : run_->second) {
    if (range.Contains(lumi)) {
      if (DEBUG) 
        PDebug("PandaLeptonicAnalyzer::PassGoodLumis",TString::Format("Accepting run=%i, lumi=%i",run,lumi));
      return true;
    }
  }

  // matched no lumi range
  if (DEBUG) 
    PDebug("PandaLeptonicAnalyzer::PassGoodLumis",TString::Format("Failing run=%i, lumi=%i",run,lumi));
  return false;
}

bool PandaLeptonicAnalyzer::PassPreselection() {

  if (preselBits==0)
    return true;

  bool isGood=false;

  if     (preselBits & kLepton) {
    if(gt->nLooseLep >= 2 && gt->looseLep1Pt > 20 && gt->looseLep2Pt > 20) isGood = true;
  }
  else if(preselBits & kLeptonFake) {
    bool passFakeTrigger = (gt->trigger & kMuFakeTrig) == kMuFakeTrig || (gt->trigger & kEGFakeTrig) == kEGFakeTrig;
    if(passFakeTrigger == true){
      double mll = 0.0;
      if(gt->nLooseLep == 2){
        TLorentzVector lep1;
        lep1.SetPtEtaPhiM(gt->looseLep1Pt,gt->looseLep1Eta,gt->looseLep1Phi,0.0);
        TLorentzVector lep2;
        lep2.SetPtEtaPhiM(gt->looseLep2Pt,gt->looseLep2Eta,gt->looseLep2Phi,0.0);
	mll = (lep1 + lep2).M();
      }
      if(mll > 70.0 || gt->nLooseLep == 1) isGood = true;
    }
  }

  return isGood;
}

void PandaLeptonicAnalyzer::CalcBJetSFs(BTagType bt, int flavor,
                double eta, double pt, double eff, double uncFactor,
                double &sf, double &sfUp, double &sfDown) 
{
  if (flavor==5) {
    sf     = btagReaders[bt]->eval_auto_bounds("central",BTagEntry::FLAV_B,eta,pt);
    sfUp   = btagReaders[bt]->eval_auto_bounds("up",BTagEntry::FLAV_B,eta,pt);
    sfDown = btagReaders[bt]->eval_auto_bounds("down",BTagEntry::FLAV_B,eta,pt);
  } else if (flavor==4) {
    sf     = btagReaders[bt]->eval_auto_bounds("central",BTagEntry::FLAV_C,eta,pt);
    sfUp   = btagReaders[bt]->eval_auto_bounds("up",BTagEntry::FLAV_C,eta,pt);
    sfDown = btagReaders[bt]->eval_auto_bounds("down",BTagEntry::FLAV_C,eta,pt);
  } else {
    sf     = btagReaders[bt]->eval_auto_bounds("central",BTagEntry::FLAV_UDSG,eta,pt);
    sfUp   = btagReaders[bt]->eval_auto_bounds("up",BTagEntry::FLAV_UDSG,eta,pt);
    sfDown = btagReaders[bt]->eval_auto_bounds("down",BTagEntry::FLAV_UDSG,eta,pt);
  }

  sfUp = uncFactor*(sfUp-sf)+sf;
  sfDown = uncFactor*(sfDown-sf)+sf;
  return;
}

void PandaLeptonicAnalyzer::EvalBTagSF(std::vector<btagcand> &cands, std::vector<double> &sfs,
               GeneralLeptonicTree::BTagShift shift,GeneralLeptonicTree::BTagJet jettype, bool do2) 
{
  float sf0 = 1, sf1 = 1, sfGT0 = 1, sf2=1;
  float prob_mc0=1, prob_data0=1;
  float prob_mc1=0, prob_data1=0;
  unsigned int nC = cands.size();

  for (unsigned int iC=0; iC!=nC; ++iC) {
    double sf_i = sfs[iC];
    double eff_i = cands[iC].eff;
    prob_mc0 *= (1-eff_i);
    prob_data0 *= (1-sf_i*eff_i);
    float tmp_mc1=1, tmp_data1=1;
    for (unsigned int jC=0; jC!=nC; ++jC) {
      if (iC==jC) continue;
      double sf_j = sfs[jC];
      double eff_j = cands[jC].eff;
      tmp_mc1 *= (1-eff_j);
      tmp_data1 *= (1-eff_j*sf_j);
    }
    prob_mc1 += eff_i * tmp_mc1;
    prob_data1 += eff_i * sf_i * tmp_data1;
  }
  
  if (nC>0) {
    sf0 = prob_data0/prob_mc0;
    sf1 = prob_data1/prob_mc1;
    sfGT0 = (1-prob_data0)/(1-prob_mc0);
  }

  GeneralLeptonicTree::BTagParams p;
  p.shift = shift;
  p.jet = jettype;
  p.tag=GeneralLeptonicTree::b0; gt->sf_btags[p] = sf0;
  p.tag=GeneralLeptonicTree::b1; gt->sf_btags[p] = sf1;
  p.tag=GeneralLeptonicTree::bGT0; gt->sf_btags[p] = sfGT0;

  if (do2) {
    float prob_mc2=0, prob_data2=0;
    unsigned int nC = cands.size();

    for (unsigned int iC=0; iC!=nC; ++iC) {
      double sf_i = sfs[iC], eff_i = cands[iC].eff;
      for (unsigned int jC=iC+1; jC!=nC; ++jC) {
        double sf_j = sfs[jC], eff_j = cands[jC].eff;
        float tmp_mc2=1, tmp_data2=1;
        for (unsigned int kC=0; kC!=nC; ++kC) {
          if (kC==iC || kC==jC) continue;
          double sf_k = sfs[kC], eff_k = cands[kC].eff;
          tmp_mc2 *= (1-eff_k);
          tmp_data2 *= (1-eff_k*sf_k);
        }
        prob_mc2 += eff_i * eff_j * tmp_mc2;
        prob_data2 += eff_i * sf_i * eff_j * sf_j * tmp_data2;
      }
    }

    if (nC>1) {
      sf2 = prob_data2/prob_mc2;
    }

    p.tag=GeneralLeptonicTree::b2; gt->sf_btags[p] = sf2;
  }

}

void PandaLeptonicAnalyzer::RegisterTrigger(TString path, std::vector<unsigned> &idxs) {
  unsigned idx = event.registerTrigger(path);
  if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::RegisterTrigger",
            TString::Format("At %u found trigger=%s",idx,path.Data()));
  idxs.push_back(idx);
}

// run
void PandaLeptonicAnalyzer::Run() {

  // INITIALIZE --------------------------------------------------------------------------
  unsigned int nEvents = tIn->GetEntries();
  unsigned int nZero = 0;
  if (lastEvent>=0 && lastEvent<(int)nEvents)
    nEvents = lastEvent;
  if (firstEvent>=0)
    nZero = firstEvent;

  if (!fOut || !tIn) {
    PError("PandaLeptonicAnalyzer::Run","NOT SETUP CORRECTLY");
    exit(1);
  }

  panda::JetCollection* jets(0);
  jets = &event.chsAK4Jets;

  // these are bins of b-tagging eff in pT and eta, derived in 8024 TT MC
  std::vector<double> vbtagpt {20.0,50.0,80.0,120.0,200.0,300.0,400.0,500.0,700.0,1000.0};
  std::vector<double> vbtageta {0.0,0.5,1.5,2.5};
  std::vector<std::vector<double>> lfeff  = {{0.081,0.065,0.060,0.063,0.072,0.085,0.104,0.127,0.162},
                       {0.116,0.097,0.092,0.099,0.112,0.138,0.166,0.185,0.222},
                       {0.173,0.145,0.149,0.175,0.195,0.225,0.229,0.233,0.250}};
  std::vector<std::vector<double>> ceff = {{0.377,0.389,0.391,0.390,0.391,0.375,0.372,0.392,0.435},
                      {0.398,0.407,0.416,0.424,0.424,0.428,0.448,0.466,0.500},
                      {0.375,0.389,0.400,0.425,0.437,0.459,0.481,0.534,0.488}};
  std::vector<std::vector<double>> beff = {{0.791,0.815,0.825,0.835,0.821,0.799,0.784,0.767,0.760},
                      {0.794,0.816,0.829,0.836,0.823,0.804,0.798,0.792,0.789},
                      {0.739,0.767,0.780,0.789,0.776,0.771,0.779,0.787,0.806}};
  Binner btagpt(vbtagpt);
  Binner btageta(vbtageta);

  JetCorrectionUncertainty *uncReader=0;
  JetCorrectionUncertainty *uncReaderAK4=0;
  FactorizedJetCorrector *scaleReaderAK4=0;

  std::vector<unsigned int> metTriggers;
  std::vector<unsigned int> phoTriggers;
  std::vector<unsigned int> muegTriggers;
  std::vector<unsigned int> mumuTriggers;
  std::vector<unsigned int> muTriggers;
  std::vector<unsigned int> muTagTriggers;
  std::vector<unsigned int> muFakeTriggers;
  std::vector<unsigned int> egegTriggers;
  std::vector<unsigned int> egTriggers;
  std::vector<unsigned int> egTagTriggers;
  std::vector<unsigned int> egFakeTriggers;

  if (1) {
    std::vector<TString> metTriggerPaths = {
          "HLT_PFMET170_NoiseCleaned",
          "HLT_PFMET170_HBHECleaned",
          "HLT_PFMET170_JetIdCleaned",
          "HLT_PFMET170_NotCleaned",
          "HLT_PFMET170_HBHE_BeamHaloCleaned",
          "HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight",
          "HLT_PFMETNoMu110_NoiseCleaned_PFMHTNoMu110_IDTight",
          "HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight",
          "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight",
          "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight",
          "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
          "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"
    };
    std::vector<TString> phoTriggerPaths = {
          "HLT_Photon175",
          "HLT_Photon165_HE10",
          "HLT_Photon36_R9Id90_HE10_IsoM",
          "HLT_Photon50_R9Id90_HE10_IsoM",
          "HLT_Photon75_R9Id90_HE10_IsoM",
          "HLT_Photon90_R9Id90_HE10_IsoM",
          "HLT_Photon120_R9Id90_HE10_IsoM",
          "HLT_Photon165_R9Id90_HE10_IsoM",
          "HLT_Photon300_NoHE",
          "HLT_ECALHT800"
    };

    std::vector<TString> muegTriggerPaths = {
	  "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
	  "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
	  "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
	  "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
	  "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
	  "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
	  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
	  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"
    };

    std::vector<TString> mumuTriggerPaths = {
	  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
	  "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
	  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
	  "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"
    };

    std::vector<TString> muTriggerPaths = {
	  "HLT_IsoMu24",
	  "HLT_IsoTkMu24",
	  "HLT_IsoMu22",
	  "HLT_IsoTkMu22",
	  "HLT_Mu45_eta2p1",
	  "HLT_Mu50"
    };

    std::vector<TString> muTagTriggerPaths = {
	  "HLT_IsoMu24",
	  "HLT_IsoTkMu24"
	  "HLT_IsoMu22",
	  "HLT_IsoTkMu22",
	  "HLT_Mu50"
    };

    std::vector<TString> muFakeTriggerPaths = {
	  "HLT_Mu8_TrkIsoVV",
	  "HLT_Mu17_TrkIsoVVL"
    };

    std::vector<TString> egegTriggerPaths = {
	  "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
	  "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf"
    };

    std::vector<TString> egTriggerPaths = {
          "HLT_Ele25_eta2p1_WPTight_Gsf",
	  "HLT_Ele27_eta2p1_WPLoose_Gsf",
	  "HLT_Ele27_WPTight_Gsf",
	  "HLT_Ele30_WPTight_Gsf",
	  "HLT_Ele35_WPLoose_Gsf",
          "HLT_Ele27_WP85_Gsf",
          "HLT_Ele27_WPLoose_Gsf",
          "HLT_Ele105_CaloIdVT_GsfTrkIdT",
          "HLT_Ele115_CaloIdVT_GsfTrkIdT",
          "HLT_Ele27_eta2p1_WPTight_Gsf",
          "HLT_Ele32_eta2p1_WPTight_Gsf",
          "HLT_ECALHT800"
    };

    std::vector<TString> egTagTriggerPaths = {
          "HLT_Ele105_CaloIdVT_GsfTrkIdT",
          "HLT_Ele115_CaloIdVT_GsfTrkIdT",
          "HLT_Ele25_eta2p1_WPTight_Gsf",
	  "HLT_Ele27_WPTight_Gsf",
          "HLT_Ele27_eta2p1_WPTight_Gsf",
	  "HLT_Ele27_eta2p1_WPLoose_Gsf"
    };

    std::vector<TString> egFakeTriggerPaths = {
          "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",
          "HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30",
          "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30"
    };

    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading MET triggers");
    for (auto path : metTriggerPaths) {
      RegisterTrigger(path,metTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading SinglePhoton triggers");
    for (auto path : phoTriggerPaths) {
      RegisterTrigger(path,phoTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading MuEG triggers");
    for (auto path : muegTriggerPaths) {
      RegisterTrigger(path,muegTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading MuMu triggers");
    for (auto path : mumuTriggerPaths) {
      RegisterTrigger(path,mumuTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading Mu triggers");
    for (auto path : muTriggerPaths) {
      RegisterTrigger(path,muTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading MuTag triggers");
    for (auto path : muTagTriggerPaths) {
      RegisterTrigger(path,muTagTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading MuFake triggers");
    for (auto path : muFakeTriggerPaths) {
      RegisterTrigger(path,muFakeTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading EGEG triggers");
    for (auto path : egegTriggerPaths) {
      RegisterTrigger(path,egegTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading EG triggers");
    for (auto path : egTriggerPaths) {
      RegisterTrigger(path,egTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading EGTag triggers");
    for (auto path : egTagTriggerPaths) {
      RegisterTrigger(path,egTagTriggers);
    }
    if (DEBUG>1) PDebug("PandaLeptonicAnalyzer::Run","Loading EGFake triggers");
    for (auto path : egFakeTriggerPaths) {
      RegisterTrigger(path,egFakeTriggers);
    }

  }

  float EGMSCALE = isData ? 1 : 1;

  // set up reporters
  unsigned int iE=0;
  ProgressReporter pr("PandaLeptonicAnalyzer::Run",&iE,&nEvents,10);
  TimeReporter tr("PandaLeptonicAnalyzer::Run",DEBUG);

  bool applyJSON = flags["applyJSON"];

  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    tr.Start();
    pr.Report();
    ResetBranches();
    event.getEntry(*tIn,iE);

    tr.TriggerEvent(TString::Format("GetEntry %u",iE));
    if (DEBUG>2) {
      PDebug("PandaLeptonicAnalyzer::Run::Dump","");
      event.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaLeptonicAnalyzer::Run::Dump","");
      event.photons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaLeptonicAnalyzer::Run::Dump","");
      event.muons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaLeptonicAnalyzer::Run::Dump","");
      event.electrons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaLeptonicAnalyzer::Run::Dump","");
      event.chsAK4Jets.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaLeptonicAnalyzer::Run::Dump","");
      event.pfMet.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("PandaLeptonicAnalyzer::Run::Dump","");
      event.metMuOnlyFix.print(std::cout, 2);
      std::cout << std::endl;
    }

    // event info
    gt->runNumber = event.runNumber;
    gt->lumiNumber = event.lumiNumber;
    gt->eventNumber = event.eventNumber;
    gt->npv = event.npv;
    gt->pu = event.npvTrue;
    gt->mcWeight = event.weight;
    gt->metFilter = (event.metFilters.pass()) ? 1 : 0;
    // these two are not need since we use muon-fixed MET
    // gt->metFilter = (gt->metFilter==1 && !event.metFilters.badMuons) ? 1 : 0;
    // gt->metFilter = (gt->metFilter==1 && !event.metFilters.duplicateMuons) ? 1 : 0;
    gt->metFilter = (gt->metFilter==1 && !event.metFilters.badPFMuons) ? 1 : 0;
    gt->metFilter = (gt->metFilter==1 && !event.metFilters.badChargedHadrons) ? 1 : 0;
    
    if(event.vertices.size() >= 1) gt->zPos = event.vertices[0].z;
    else                           gt->zPos = 0.0;

    // save triggers
    for (auto iT : metTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kMETTrig;
      break;
     }
    }
    for (auto iT : phoTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kSinglePhoTrig;
      break;
     }
    }
    for (auto iT : muegTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kMuEGTrig;
      break;
     }
    }
    for (auto iT : mumuTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kMuMuTrig;
      break;
     }
    }
    for (auto iT : muTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kMuTrig;
      break;
     }
    }
    for (auto iT : muTagTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kMuTagTrig;
      break;
     }
    }
    for (auto iT : muFakeTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kMuFakeTrig;
      break;
     }
    }
    for (auto iT : egegTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kEGEGTrig;
      break;
     }
    }
    for (auto iT : egTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kEGTrig;
      break;
     }
    }
    for (auto iT : egTagTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kEGTagTrig;
      break;
     }
    }
    for (auto iT : egFakeTriggers) {
     if (event.triggerFired(iT)) {
      gt->trigger |= kEGFakeTrig;
      break;
     }
    }

    if (isData) {
      // check the json
      if (applyJSON && !PassGoodLumis(gt->runNumber,gt->lumiNumber))
        continue;

    } else {
      gt->sf_pu     = GetCorr(cPU    ,gt->pu);
      gt->sf_puUp   = GetCorr(cPUUp  ,gt->pu);
      gt->sf_puDown = GetCorr(cPUDown,gt->pu);
    }

    if (uncReader==0) {
      if (isData) {
        TString thisEra = eras.getEra(gt->runNumber);
        for (auto &iter : ak8UncReader) {
          if (! iter.first.Contains("data"))
            continue;
          if (iter.first.Contains(thisEra)) {
            uncReader = iter.second;
            uncReaderAK4 = ak4UncReader[iter.first];
            scaleReaderAK4 = ak4ScaleReader[iter.first];
            break;
          }
        }
      } else {
        uncReader = ak8UncReader["MC"];
        uncReaderAK4 = ak4UncReader["MC"];
        scaleReaderAK4 = ak4ScaleReader["MC"];
      }
    }

    tr.TriggerEvent("initialize");

    // met
    gt->pfmet = event.pfMet.pt;
    gt->pfmetphi = event.pfMet.phi;
    gt->pfmetRaw = event.rawMet.pt;
    gt->pfmetUp = event.pfMet.ptCorrUp;
    gt->pfmetDown = event.pfMet.ptCorrDown;
    gt->puppimet = event.puppiMet.pt;
    gt->puppimetphi = event.puppiMet.phi;
    gt->calomet = event.caloMet.pt;
    gt->calometphi = event.caloMet.phi;
    gt->trkmet = event.trkMet.pt;
    gt->trkmetphi = event.trkMet.phi;
    TLorentzVector vPFMET, vPuppiMET;
    vPFMET.SetPtEtaPhiM(gt->pfmet,0,gt->pfmetphi,0);
    vPuppiMET.SetPtEtaPhiM(gt->puppimet,0,gt->puppimetphi,0);
    TVector2 vMETNoMu; vMETNoMu.SetMagPhi(gt->pfmet,gt->pfmetphi); //       for trigger eff

    tr.TriggerEvent("met");

    //electrons
    std::vector<panda::Lepton*> looseLeps;
    for (auto& ele : event.electrons) {
     float pt = ele.pt()*EGMSCALE; float eta = ele.eta(); float aeta = fabs(eta);
      if (pt<10 || aeta>2.5)
        continue;
      if (!ele.veto)
        continue;
      looseLeps.push_back(&ele);
      matchLeps.push_back(&ele);
    }

    // muons
    for (auto& mu : event.muons) {
     float pt = mu.pt(); float eta = mu.eta(); float aeta = fabs(eta);
      if (pt<10 || aeta>2.4)
        continue;
      if (!mu.loose)
        continue;
      looseLeps.push_back(&mu);
      matchLeps.push_back(&mu);
      TVector2 vMu; vMu.SetMagPhi(pt,mu.phi());
      vMETNoMu += vMu;
    }
    gt->pfmetnomu = vMETNoMu.Mod();

    // now consider all leptons
    gt->nLooseLep = looseLeps.size();
    if (gt->nLooseLep>0) {
     auto ptsort([](panda::Lepton const* l1, panda::Lepton const* l2)->bool {
       return l1->pt() > l2->pt();
      });
     int nToSort = gt->nLooseLep;
     std::partial_sort(looseLeps.begin(),looseLeps.begin()+nToSort,looseLeps.end(),ptsort);
    }
    int lep_counter=1;
    for (auto* lep : looseLeps) {
      if      (lep_counter==1) {
       gt->looseLep1Pt  = lep->pt();
       gt->looseLep1Eta = lep->eta();
       gt->looseLep1Phi = lep->phi();
      }
      else if (lep_counter==2) {
       gt->looseLep2Pt  = lep->pt();
       gt->looseLep2Eta = lep->eta();
       gt->looseLep2Phi = lep->phi();
      } 
      else if (lep_counter==3) {
       gt->looseLep3Pt  = lep->pt();
       gt->looseLep3Eta = lep->eta();
       gt->looseLep3Phi = lep->phi();
      } 
      else if (lep_counter==4) {
       gt->looseLep4Pt  = lep->pt();
       gt->looseLep4Eta = lep->eta();
       gt->looseLep4Phi = lep->phi();
      } 
      else {
        break;
      }
      // now specialize lepton types
      panda::Muon *mu = dynamic_cast<panda::Muon*>(lep);
      if (mu!=NULL) {
        bool isLoose  = mu->loose;
        bool isFake   = mu->tight  && mu->combIso()/mu->pt() < 0.4 && mu->chIso/mu->pt() < 0.4;
        bool isMedium = (mu->medium || mu-> mediumBtoF) && mu->combIso()/mu->pt() < 0.15;
        bool isTight  = mu->tight  && mu->combIso()/mu->pt() < 0.15;
        bool isDxyz   = MuonIP(mu->dxy,mu->dz);
	bool isTrigger = false;

        for (auto& trigObj : event.triggerObjects.filterObjects("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09")) { // HLT_IsoMu24
         if (trigObj->pt() > 0 && mu->pt() > 0 && trigObj->dR(*mu) < 0.1) isTrigger = true;
        }
	if(isTrigger == false){
          for (auto& trigObj : event.triggerObjects.filterObjects("hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09")) { // HLT_IsoTkMu24
           if (trigObj->pt() > 0 && mu->pt() > 0 && trigObj->dR(*mu) < 0.1) isTrigger = true;
          }
	}
	if(isTrigger == false){
          for (auto& trigObj : event.triggerObjects.filterObjects("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09")) { // HLT_IsoMu22
           if (trigObj->pt() > 0 && mu->pt() > 0 && trigObj->dR(*mu) < 0.1) isTrigger = true;
          }
	}
	if(isTrigger == false){
          for (auto& trigObj : event.triggerObjects.filterObjects("hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09")) { // HLT_IsoTkMu22
           if (trigObj->pt() > 0 && mu->pt() > 0 && trigObj->dR(*mu) < 0.1) isTrigger = true;
          }
	}
	if(isTrigger == false){
          for (auto& trigObj : event.triggerObjects.filterObjects("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q")) { // HLT_Mu50
           if (trigObj->pt() > 0 && mu->pt() > 0 && trigObj->dR(*mu) < 0.1) isTrigger = true;
          }
	}

        isTrigger = isTrigger || mu->triggerMatch[panda::Muon::fIsoMu24] || mu->triggerMatch[panda::Muon::fIsoTkMu24] || mu->triggerMatch[panda::Muon::fIsoMu22er] || mu->triggerMatch[panda::Muon::fIsoTkMu22er] || mu->triggerMatch[panda::Muon::fMu50];

	bool isFakeTrigger = false;
        for (auto& trigObj : event.triggerObjects.filterObjects("hltL3fL1sMu5L1f0L2f5L3Filtered8TkIsoFiltered0p4")) { // HLT_Mu8_TrkIsoVVL
         if (trigObj->pt() > 0 && mu->pt() > 0 && trigObj->dR(*mu) < 0.1) isFakeTrigger = true;
        }
	if(isFakeTrigger == false){
          for (auto& trigObj : event.triggerObjects.filterObjects("hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4")) { // HLT_Mu17_TrkIsoVVL
           if (trigObj->pt() > 0 && mu->pt() > 0 && trigObj->dR(*mu) < 0.1) isFakeTrigger = true;
          }
	}

        if      (lep_counter==1) {
          gt->looseLep1SCEta = mu->pfPt; // sure, it is a hack!
          gt->looseLep1RegPt = mu->trkLayersWithMmt; // sure, it is a hack!
          gt->looseLep1SmePt = mu->pt();
          gt->looseLep1PdgId = mu->charge*-13;
          if(isLoose)      gt->looseLep1SelBit |= kLoose;
          if(isFake)       gt->looseLep1SelBit |= kFake;
          if(isMedium)     gt->looseLep1SelBit |= kMedium;
          if(isTight)      gt->looseLep1SelBit |= kTight;
          if(isDxyz)       gt->looseLep1SelBit |= kDxyz;
          if(isTrigger)    gt->looseLep1SelBit |= kTrigger;
          if(isFakeTrigger)gt->looseLep1SelBit |= kFakeTrigger;
        }
	else if (lep_counter==2) {
          gt->looseLep2SCEta = mu->pfPt; // sure, it is a hack!
          gt->looseLep2RegPt = mu->trkLayersWithMmt; // sure, it is a hack!
          gt->looseLep2SmePt = mu->pt();
          gt->looseLep2PdgId = mu->charge*-13;
          if(isLoose)      gt->looseLep2SelBit |= kLoose;
          if(isFake)       gt->looseLep2SelBit |= kFake;
          if(isMedium)     gt->looseLep2SelBit |= kMedium;
          if(isTight)      gt->looseLep2SelBit |= kTight;
          if(isDxyz)       gt->looseLep2SelBit |= kDxyz;
          if(isTrigger)    gt->looseLep2SelBit |= kTrigger;
          if(isFakeTrigger)gt->looseLep2SelBit |= kFakeTrigger;
        }
	else if (lep_counter==3) {
          gt->looseLep3SCEta = mu->pfPt; // sure, it is a hack!
          gt->looseLep3RegPt = mu->trkLayersWithMmt; // sure, it is a hack!
          gt->looseLep3SmePt = mu->pt();
          gt->looseLep3PdgId = mu->charge*-13;
          if(isLoose)      gt->looseLep3SelBit |= kLoose;
          if(isFake)       gt->looseLep3SelBit |= kFake;
          if(isMedium)     gt->looseLep3SelBit |= kMedium;
          if(isTight)      gt->looseLep3SelBit |= kTight;
          if(isDxyz)       gt->looseLep3SelBit |= kDxyz;
          if(isTrigger)    gt->looseLep3SelBit |= kTrigger;
          if(isFakeTrigger)gt->looseLep3SelBit |= kFakeTrigger;
        }
	else if (lep_counter==4) {
          gt->looseLep4SCEta = mu->pfPt; // sure, it is a hack!
          gt->looseLep4RegPt = mu->trkLayersWithMmt; // sure, it is a hack!
          gt->looseLep4SmePt = mu->pt();
          gt->looseLep4PdgId = mu->charge*-13;
          if(isLoose)      gt->looseLep4SelBit |= kLoose;
          if(isFake)       gt->looseLep4SelBit |= kFake;
          if(isMedium)     gt->looseLep4SelBit |= kMedium;
          if(isTight)      gt->looseLep4SelBit |= kTight;
          if(isDxyz)       gt->looseLep4SelBit |= kDxyz;
          if(isTrigger)    gt->looseLep4SelBit |= kTrigger;
          if(isFakeTrigger)gt->looseLep4SelBit |= kFakeTrigger;
        }
      } else {
        panda::Electron *ele = dynamic_cast<panda::Electron*>(lep);
        bool isLoose  = ele->loose;
        bool isFake   = ele->hltsafe;
        bool isMedium = ele->medium;
        bool isTight  = ele->tight;
        bool isDxyz   = ElectronIP(ele->eta(),ele->dxy,ele->dz);
	bool isTrigger = false;

        for (auto& trigObj : event.triggerObjects.filterObjects("hltEle105CaloIdVTGsfTrkIdTGsfDphiFilter")) { // HLT_Ele105_CaloIdVT_GsfTrkIdT
         if (trigObj->pt() > 0 && ele->pt() > 0 && trigObj->dR(*ele) < 0.1) isTrigger = true;
        }
	if(isTrigger == false){
          for (auto& trigObj : event.triggerObjects.filterObjects("hltEle115CaloIdVTGsfTrkIdTGsfDphiFilter")) { // HLT_Ele115_CaloIdVT_GsfTrkIdT
           if (trigObj->pt() > 0 && ele->pt() > 0 && trigObj->dR(*ele) < 0.1) isTrigger = true;
          }
	}
	if(isTrigger == false){
          for (auto& trigObj : event.triggerObjects.filterObjects("hltEle25erWPTightGsfTrackIsoFilter")) { // HLT_Ele25_eta2p1_WPTight_Gsf
           if (trigObj->pt() > 0 && ele->pt() > 0 && trigObj->dR(*ele) < 0.1) isTrigger = true;
          }
	}
	if(isTrigger == false){
          for (auto& trigObj : event.triggerObjects.filterObjects("hltEle27WPTightGsfTrackIsoFilter")) { // HLT_Ele27_WPTight_Gsf
           if (trigObj->pt() > 0 && ele->pt() > 0 && trigObj->dR(*ele) < 0.1) isTrigger = true;
          }
	}
	if(isTrigger == false){
          for (auto& trigObj : event.triggerObjects.filterObjects("hltEle27erWPTightGsfTrackIsoFilter")) { // HLT_Ele27_eta2p1_WPTight_Gsf
           if (trigObj->pt() > 0 && ele->pt() > 0 && trigObj->dR(*ele) < 0.1) isTrigger = true;
          }
	}
	if(isTrigger == false){
          for (auto& trigObj : event.triggerObjects.filterObjects("hltEle27erWPLooseGsfTrackIsoFilter")) { // HLT_Ele27_eta2p1_WPLoose_Gsf
           if (trigObj->pt() > 0 && ele->pt() > 0 && trigObj->dR(*ele) < 0.1) isTrigger = true;
          }
	}

        isTrigger = isTrigger || ele->triggerMatch[panda::Electron::fEl25Tight] || ele->triggerMatch[panda::Electron::fEl27Tight] || ele->triggerMatch[panda::Electron::fEl27Loose];

	bool isFakeTrigger = false;
        for (auto& trigObj : event.triggerObjects.filterObjects("hltEle12PFJet30EleCleaned")) { // HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30
         if (trigObj->pt() > 0 && ele->pt() > 0 && trigObj->dR(*ele) < 0.1) isFakeTrigger = true;
        }
	if(isFakeTrigger == false){
          for (auto& trigObj : event.triggerObjects.filterObjects("hltEle17PFJet30EleCleaned")) { // HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30
           if (trigObj->pt() > 0 && ele->pt() > 0 && trigObj->dR(*ele) < 0.1) isFakeTrigger = true;
          }
	}
	if(isFakeTrigger == false){
          for (auto& trigObj : event.triggerObjects.filterObjects("hltEle23PFJet30EleCleaned")) { // HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30
           if (trigObj->pt() > 0 && ele->pt() > 0 && trigObj->dR(*ele) < 0.1) isFakeTrigger = true;
          }
	}

	if(TMath::Abs(ele->eta()-ele->superCluster->eta) > 0.2) printf("Potential issue ele/sc: dist: %f - %f/%f/%f vs. %f/%f/%f\n",TMath::Abs(ele->eta()-ele->superCluster->eta),ele->pt(),ele->eta(),ele->phi(),ele->superCluster->rawPt,ele->superCluster->eta,ele->superCluster->phi);
        if      (lep_counter==1) {
          gt->looseLep1SCEta = ele->superCluster->eta;
          gt->looseLep1RegPt = ele->regPt;
          gt->looseLep1SmePt = ele->smearedPt;
          gt->looseLep1Pt *= EGMSCALE;
          gt->looseLep1PdgId = ele->charge*-11;
          if(isLoose)      gt->looseLep1SelBit |= kLoose;
          if(isFake)       gt->looseLep1SelBit |= kFake;
          if(isMedium)     gt->looseLep1SelBit |= kMedium;
          if(isTight)      gt->looseLep1SelBit |= kTight;
          if(isDxyz)       gt->looseLep1SelBit |= kDxyz;
          if(isTrigger)    gt->looseLep1SelBit |= kTrigger;
          if(isFakeTrigger)gt->looseLep1SelBit |= kFakeTrigger;
        }
	else if (lep_counter==2) {
          gt->looseLep2SCEta = ele->superCluster->eta;
          gt->looseLep2RegPt = ele->regPt;
          gt->looseLep2SmePt = ele->smearedPt;
          gt->looseLep2Pt *= EGMSCALE;
          gt->looseLep2PdgId = ele->charge*-11;
          if(isLoose)      gt->looseLep2SelBit |= kLoose;
          if(isFake)       gt->looseLep2SelBit |= kFake;
          if(isMedium)     gt->looseLep2SelBit |= kMedium;
          if(isTight)      gt->looseLep2SelBit |= kTight;
          if(isDxyz)       gt->looseLep2SelBit |= kDxyz;
          if(isTrigger)    gt->looseLep2SelBit |= kTrigger;
          if(isFakeTrigger)gt->looseLep2SelBit |= kFakeTrigger;
        }
	else if (lep_counter==3) {
          gt->looseLep3SCEta = ele->superCluster->eta;
          gt->looseLep3RegPt = ele->regPt;
          gt->looseLep3SmePt = ele->smearedPt;
          gt->looseLep3Pt *= EGMSCALE;
          gt->looseLep3PdgId = ele->charge*-11;
          if(isLoose)      gt->looseLep3SelBit |= kLoose;
          if(isFake)       gt->looseLep3SelBit |= kFake;
          if(isMedium)     gt->looseLep3SelBit |= kMedium;
          if(isTight)      gt->looseLep3SelBit |= kTight;
          if(isDxyz)       gt->looseLep3SelBit |= kDxyz;
          if(isTrigger)    gt->looseLep3SelBit |= kTrigger;
          if(isFakeTrigger)gt->looseLep3SelBit |= kFakeTrigger;
        }
	else if (lep_counter==4) {
          gt->looseLep4SCEta = ele->superCluster->eta;
          gt->looseLep4RegPt = ele->regPt;
          gt->looseLep4SmePt = ele->smearedPt;
          gt->looseLep4Pt *= EGMSCALE;
          gt->looseLep4PdgId = ele->charge*-11;
          if(isLoose)      gt->looseLep4SelBit |= kLoose;
          if(isFake)       gt->looseLep4SelBit |= kFake;
          if(isMedium)     gt->looseLep4SelBit |= kMedium;
          if(isTight)      gt->looseLep4SelBit |= kTight;
          if(isDxyz)       gt->looseLep4SelBit |= kDxyz;
          if(isTrigger)    gt->looseLep4SelBit |= kTrigger;
          if(isFakeTrigger)gt->looseLep4SelBit |= kFakeTrigger;
        }
      }
      ++lep_counter;
    }

    tr.TriggerEvent("leptons");

    // photons
    gt->nLoosePhoton = 0;
    for (auto& pho : event.photons) {
      if (!pho.medium || !pho.pixelVeto)
        continue;
      float pt = pho.pt() * EGMSCALE;
      if (pt<1) continue;
      float eta = pho.eta(), phi = pho.phi();
      if (pt<20 || fabs(eta)>2.5)
        continue;
      if (IsMatched(&matchLeps,0.16,eta,phi))
        continue;
      gt->nLoosePhoton++;
      matchPhos.push_back(&pho);
      if (gt->nLoosePhoton==1) {
        gt->loosePho1Pt = pt;
        gt->loosePho1Eta = eta;
        gt->loosePho1Phi = phi;
      }
    }

    tr.TriggerEvent("photons");

    // first identify interesting jets
    vector<panda::Jet*> cleaned30Jets,cleaned20Jets;
    vector<int> btagindices;
    TLorentzVector vJet;
    panda::Jet *jet1=0, *jet2=0, *jet3=0, *jet4=0;
    panda::Jet *jetUp1=0, *jetUp2=0, *jetUp3=0, *jetUp4=0;
    panda::Jet *jetDown1=0, *jetDown2=0, *jetDown3=0, *jetDown4=0;
    gt->dphipuppimet=999; gt->dphipfmet=999;
    float maxJetEta = 4.7;
    unsigned nJetDPhi = 1;

    for (auto& jet : *jets) {

      // only do eta-phi checks here
      if (abs(jet.eta()) > maxJetEta)
         continue;
      if (IsMatched(&matchLeps,0.16,jet.eta(),jet.phi()))
         continue;

      bool isLoose = jet.loose;
      bool isTight = jet.tight;

      if (jet.pt()>20 && jet.csv>0.5426) ++(gt->jetNLBtags);
      if (jet.pt()>20 && jet.csv>0.8484) ++(gt->jetNMBtags);
      if (jet.pt()>20 && jet.csv>0.9535) ++(gt->jetNTBtags);

      if (jet.pt()>20) cleaned20Jets.push_back(&jet); // to be used for btagging SFs

      if (jet.pt()>30) { // nominal jets
	cleaned30Jets.push_back(&jet);
	if      (cleaned30Jets.size()==1) {
          jet1 = &jet;
          gt->jet1Pt   = jet.pt();
          gt->jet1Eta  = jet.eta();
          gt->jet1Phi  = jet.phi();
          gt->jet1BTag = jet.csv;
          if(isLoose) gt->jet1SelBit |= kLoose;
          if(isTight) gt->jet1SelBit |= kTight;
	}
	else if (cleaned30Jets.size()==2) {
          jet2 = &jet;
          gt->jet2Pt   = jet.pt();
          gt->jet2Eta  = jet.eta();
          gt->jet2Phi  = jet.phi();
          gt->jet2BTag = jet.csv;
          if(isLoose) gt->jet2SelBit |= kLoose;
          if(isTight) gt->jet2SelBit |= kTight;
	}
	else if (cleaned30Jets.size()==3) {
          jet3 = &jet;
          gt->jet3Pt   = jet.pt();
          gt->jet3Eta  = jet.eta();
          gt->jet3Phi  = jet.phi();
          gt->jet3BTag = jet.csv;
          if(isLoose) gt->jet3SelBit |= kLoose;
          if(isTight) gt->jet3SelBit |= kTight;
	}
	else if (cleaned30Jets.size()==4) {
          jet4 = &jet;
          gt->jet4Pt   = jet.pt();
          gt->jet4Eta  = jet.eta();
          gt->jet4Phi  = jet.phi();
          gt->jet4BTag = jet.csv;
          if(isLoose) gt->jet4SelBit |= kLoose;
          if(isTight) gt->jet4SelBit |= kTight;
	}

	// compute dphi wrt mets
	if (cleaned30Jets.size() <= nJetDPhi) {
          vJet.SetPtEtaPhiM(jet.pt(),jet.eta(),jet.phi(),jet.m());
          gt->dphipuppimet = std::min(fabs(vJet.DeltaPhi(vPuppiMET)),(double)gt->dphipuppimet);
          gt->dphipfmet = std::min(fabs(vJet.DeltaPhi(vPFMET)),(double)gt->dphipfmet);
	}
      }

      // do jes variation OUTSIDE of pt>30 check
      if (jet.ptCorrUp>30) {
	if (jet.ptCorrUp > gt->jet1PtUp) {
          if (jetUp1) {
            jetUp4 = jetUp3;
            gt->jet4PtUp  = gt->jet3PtUp;
            gt->jet4EtaUp = gt->jet3EtaUp;
            jetUp3 = jetUp2;
            gt->jet3PtUp  = gt->jet2PtUp;
            gt->jet3EtaUp = gt->jet2EtaUp;
            jetUp2 = jetUp1;
            gt->jet2PtUp  = gt->jet1PtUp;
            gt->jet2EtaUp = gt->jet1EtaUp;
          }
          jetUp1 = &jet;
          gt->jet1PtUp  = jet.ptCorrUp;
          gt->jet1EtaUp = jet.eta();
	}
	else if (jet.ptCorrUp > gt->jet2PtUp) {
          if (jetUp2) {
            jetUp4 = jetUp3;
            gt->jet4PtUp  = gt->jet3PtUp;
            gt->jet4EtaUp = gt->jet3EtaUp;
            jetUp3 = jetUp2;
            gt->jet3PtUp  = gt->jet2PtUp;
            gt->jet3EtaUp = gt->jet2EtaUp;
          }
          jetUp2 = &jet;
          gt->jet2PtUp  = jet.ptCorrUp;
          gt->jet2EtaUp = jet.eta();
	}
	else if (jet.ptCorrUp > gt->jet3PtUp) {
          if (jetUp3) {
            jetUp4 = jetUp3;
            gt->jet4PtUp  = gt->jet3PtUp;
            gt->jet4EtaUp = gt->jet3EtaUp;
          }
          jetUp3 = &jet;
          gt->jet3PtUp  = jet.ptCorrUp;
          gt->jet3EtaUp = jet.eta();
	}
	else if (jet.ptCorrUp > gt->jet4PtUp) {
          jetUp4 = &jet;
          gt->jet4PtUp  = jet.ptCorrUp;
          gt->jet4EtaUp = jet.eta();
	}
      }

      if (jet.ptCorrDown>30) {
	if (jet.ptCorrDown > gt->jet1PtDown) {
          if (jetDown1) {
            jetDown4 = jetDown3;
            gt->jet4PtDown  = gt->jet3PtDown;
            gt->jet4EtaDown = gt->jet3EtaDown;
            jetDown3 = jetDown2;
            gt->jet3PtDown  = gt->jet2PtDown;
            gt->jet3EtaDown = gt->jet2EtaDown;
            jetDown2 = jetDown1;
            gt->jet2PtDown  = gt->jet1PtDown;
            gt->jet2EtaDown = gt->jet1EtaDown;
          }
          jetDown1 = &jet;
          gt->jet1PtDown  = jet.ptCorrDown;
          gt->jet1EtaDown = jet.eta();
	}
	else if (jet.ptCorrDown > gt->jet2PtDown) {
          if (jetDown2) {
            jetDown4 = jetDown3;
            gt->jet4PtDown  = gt->jet3PtDown;
            gt->jet4EtaDown = gt->jet3EtaDown;
            jetDown3 = jetDown2;
            gt->jet3PtDown  = gt->jet2PtDown;
            gt->jet3EtaDown = gt->jet2EtaDown;
          }
          jetDown2 = &jet;
          gt->jet2PtDown  = jet.ptCorrDown;
          gt->jet2EtaDown = jet.eta();
	}
	else if (jet.ptCorrDown > gt->jet3PtDown) {
          if (jetDown3) {
            jetDown4 = jetDown3;
            gt->jet4PtDown  = gt->jet3PtDown;
            gt->jet4EtaDown = gt->jet3EtaDown;
          }
          jetDown3 = &jet;
          gt->jet3PtDown  = jet.ptCorrDown;
          gt->jet3EtaDown = jet.eta();
	}
	else if (jet.ptCorrDown > gt->jet4PtDown) {
          jetDown4 = &jet;
          gt->jet4PtDown  = jet.ptCorrDown;
          gt->jet4EtaDown = jet.eta();
	}
      }
    } // Jet loop

    gt->nJet = cleaned30Jets.size();

    tr.TriggerEvent("jets");

    for (auto& tau : event.taus) {
      if (!tau.decayMode || !tau.decayModeNew)
        continue;
      if (!tau.looseIsoMVA)
        continue;
      if (tau.pt()<18 || fabs(tau.eta())>2.3)
        continue;
      if (IsMatched(&matchLeps,0.16,tau.eta(),tau.phi()))
        continue;
      gt->nTau++;
    }

    tr.TriggerEvent("taus");

    // scale and PDF weights, if they exist
    gt->pdfUp = 1; gt->pdfDown = 1;
    if (!isData) {
      gt->pdfUp = 1 + event.genReweight.pdfDW;
      gt->pdfDown = 1 - event.genReweight.pdfDW;
      auto &genReweight = event.genReweight;
      for (unsigned iS=0; iS!=6; ++iS) {
        float s=1;
        switch (iS) {
          case 0:
            s = genReweight.r1f2DW; break;
          case 1:
            s = genReweight.r1f5DW; break;
          case 2:
            s = genReweight.r2f1DW; break;
          case 3:
            s = genReweight.r2f2DW; break;
          case 4:
            s = genReweight.r5f1DW; break;
          case 5:
            s = genReweight.r5f5DW; break;
          default:
            break;
        }
        gt->scale[iS] = s; 
      }
      tr.TriggerEvent("qcd uncertainties");

      unsigned nW = wIDs.size();
      if (nW) {
        for (unsigned iW=0; iW!=nW; ++iW) {
          gt->signal_weights[wIDs[iW]] = event.genReweight.genParam[iW];
        }
      }
    }

    double maxQCDscale = 1.0;
    if(gt->scale[0] != -1){
      maxQCDscale = (TMath::Abs(1+gt->scale[0])+TMath::Abs(1+gt->scale[1])+TMath::Abs(1+gt->scale[2])+
    		     TMath::Abs(1+gt->scale[3])+TMath::Abs(1+gt->scale[4])+TMath::Abs(1+gt->scale[4]))/6.0;
    }

    gt->sf_tt = 1;
    gt->genLep1Pt = 0;
    gt->genLep1Eta = -1;
    gt->genLep1Phi = -1;
    gt->genLep1PdgId = 0;
    gt->genLep2Pt = 0;
    gt->genLep2Eta = -1;
    gt->genLep2Phi = -1;
    gt->genLep2PdgId = 0;
    gt->looseGenLep1PdgId = 0;
    gt->looseGenLep2PdgId = 0;
    gt->looseGenLep3PdgId = 0;
    gt->looseGenLep4PdgId = 0;
    TLorentzVector v1,v2,v3,v4;
    if (gt->nLooseLep>=1) {
      panda::Lepton *lep1=looseLeps[0];
      v1.SetPtEtaPhiM(lep1->pt(),lep1->eta(),lep1->phi(),lep1->m());
    }
    if (gt->nLooseLep>=2) {
      panda::Lepton *lep2=looseLeps[1];
      v2.SetPtEtaPhiM(lep2->pt(),lep2->eta(),lep2->phi(),lep2->m());
    }
    if (gt->nLooseLep>=3) {
      panda::Lepton *lep3=looseLeps[2];
      v3.SetPtEtaPhiM(lep3->pt(),lep3->eta(),lep3->phi(),lep3->m());
    }
    if (gt->nLooseLep>=4) {
      panda::Lepton *lep4=looseLeps[3];
      v4.SetPtEtaPhiM(lep4->pt(),lep4->eta(),lep4->phi(),lep4->m());
    }
    // gen lepton matching
    if (!isData) {
      TLorentzVector thePartonZ(0,0,0,0); int nPartonLeptons = 0;
      int nPartons = event.partons.size();
      for (int iG=0; iG!=nPartons; ++iG) {
        auto& part(event.partons.at(iG));
        unsigned int abspdgid = abs(part.pdgid);
        if (abspdgid == 11 || abspdgid == 13 || abspdgid == 15) {
          nPartonLeptons++;
          TLorentzVector partonLepton;
          partonLepton.SetPtEtaPhiM(part.pt(),part.eta(),part.phi(),part.m());
          thePartonZ = thePartonZ + partonLepton;        
	}
      }
      //printf("kZPtCut %d %f\n",nPartonLeptons,thePartonZ.Pt());
      if(processType == kZPtCut && nPartonLeptons == 2 && thePartonZ.Pt() > 100) continue;

      std::vector<int> targetsLepton;
      std::vector<int> targetsPhoton;
      std::vector<int> targetsV;
      std::vector<int> targetsTop;
      std::vector<int> targetsN;
      int nGen = event.genParticles.size();
      for (int iG=0; iG!=nGen; ++iG) {
        auto& part(event.genParticles.at(iG));
        int pdgid = part.pdgid;
        unsigned int abspdgid = abs(pdgid);
        if ((abspdgid == 11 || abspdgid == 13) && part.finalState &&
	    (part.testFlag(GenParticle::kIsPrompt) || part.statusFlags == GenParticle::kIsPrompt ||
	     part.testFlag(GenParticle::kIsTauDecayProduct) || part.testFlag(GenParticle::kIsPromptTauDecayProduct) || 
	     part.testFlag(GenParticle::kIsDirectTauDecayProduct) || part.testFlag(GenParticle::kIsDirectPromptTauDecayProduct) ||
	     (part.parent.isValid() && abs(part.parent->pdgid) == 15)))
          targetsLepton.push_back(iG);

        if (abspdgid == 22 && part.finalState)
          targetsPhoton.push_back(iG);

        if (abspdgid == 23 || abspdgid == 24)
          targetsV.push_back(iG);

        if (abspdgid == 6)
          targetsTop.push_back(iG);

        if (abspdgid == 12 || abspdgid == 14 || abspdgid == 16)
          targetsN.push_back(iG);

      } //looking for targets

      std::vector<int> targetsJet;
      int nGoodGenJets[3] = {0,0,0};
      for (int iG=0; iG!=event.ak4GenJets.size(); ++iG) {
        auto& partj(event.ak4GenJets.at(iG));
        if(TMath::Abs(partj.eta()) >= 5.0) continue;
	if(partj.pt() <= 25) continue;

        bool isLepton = kFALSE;
        for (int iG : targetsLepton) {
          auto& partl(event.genParticles.at(iG));
          if(TMath::Abs(partl.eta()) >= 2.5) continue;
	  if(partl.pt() <= 20) continue;

//           // check there is no further copy:
//           bool isLastCopy=true;
//           for (int kG : targetsLepton) {
//             if (event.genParticles.at(kG).parent.isValid() && event.genParticles.at(kG).parent.get() == &partl) {
//               isLastCopy=false;
//               break;
//             }
//           }
//           if (!isLastCopy)
//             continue;

	  if(DeltaR2(partj.eta(),partj.phi(),partl.eta(),partl.phi()) < 0.09){
	     isLepton = kTRUE;
          }
	}
	if(isLepton == kTRUE) continue;
        if(partj.pt() > 25) {nGoodGenJets[0]++;}
        if(partj.pt() > 30) {nGoodGenJets[1]++; targetsJet.push_back(iG);}
        if(partj.pt() > 35) {nGoodGenJets[2]++;}
      }

      double deltaEtaJJGen = 0; double mJJGen = 0.;
      if(targetsJet.size() >= 2) {
        auto& partj1(event.ak4GenJets.at(0));
        auto& partj2(event.ak4GenJets.at(1));
        TLorentzVector genjet1;
        genjet1.SetPtEtaPhiM(partj1.pt(),partj1.eta(),partj1.phi(),0.0);
        TLorentzVector genjet2;
        genjet2.SetPtEtaPhiM(partj2.pt(),partj2.eta(),partj2.phi(),0.0);
	TLorentzVector dijet = genjet1 + genjet2;
	deltaEtaJJGen = TMath::Abs(partj1.eta()-partj2.eta());
	mJJGen = dijet.M();
      }

      TLorentzVector the_rhoP4(0,0,0,0);
      double bosonPtMin = 1000000000;
      for (int iG : targetsLepton) {
        auto& part(event.genParticles.at(iG));
        TLorentzVector dressedLepton;
        dressedLepton.SetPtEtaPhiM(part.pt(),part.eta(),part.phi(),part.m());

//         // check there is no further copy:
//         bool isLastCopy=true;
//         for (int kG : targetsLepton) {
//           if (event.genParticles.at(kG).parent.isValid() && event.genParticles.at(kG).parent.get() == &part) {
//             isLastCopy=false;
//             break;
//           }
//         }
//         if (!isLastCopy)
//           continue;
	
	the_rhoP4 = the_rhoP4 + dressedLepton;
	
        for (int jG : targetsPhoton) {
          auto& partj(event.genParticles.at(jG));

//           // check there is no further copy:
//           bool isLastCopy=true;
//           for (int kG : targetsPhoton) {
//             if (event.genParticles.at(kG).parent.isValid() && event.genParticles.at(kG).parent.get() == &part) {
//               isLastCopy=false;
//               break;
//             }
//           }
//           if (!isLastCopy)
//             continue;

	  if(abs(partj.pdgid) == 22 && DeltaR2(part.eta(),part.phi(),partj.eta(),partj.phi()) < 0.01) {
            TLorentzVector photonV;
            photonV.SetPtEtaPhiM(partj.pt(),partj.eta(),partj.phi(),partj.m());
	    dressedLepton += photonV;
	  }
	}

        int thePdgId = part.pdgid;
	if (part.testFlag(GenParticle::kIsTauDecayProduct) || part.testFlag(GenParticle::kIsPromptTauDecayProduct) || 
	    part.testFlag(GenParticle::kIsDirectTauDecayProduct) || part.testFlag(GenParticle::kIsDirectPromptTauDecayProduct) ||
	    (part.parent.isValid() && abs(part.parent->pdgid) == 15)) thePdgId = 15 * part.pdgid/abs(part.pdgid);
	if     (dressedLepton.Pt() > gt->genLep1Pt){
          gt->genLep2Pt    = gt->genLep1Pt; 
          gt->genLep2Eta   = gt->genLep1Eta;
          gt->genLep2Phi   = gt->genLep1Phi;
          gt->genLep2PdgId = gt->genLep1PdgId; 
          gt->genLep1Pt    = dressedLepton.Pt();
          gt->genLep1Eta   = dressedLepton.Eta();
          gt->genLep1Phi   = dressedLepton.Phi();
          gt->genLep1PdgId = thePdgId;
        }
	else if(dressedLepton.Pt() > gt->genLep2Pt){
          gt->genLep2Pt    = dressedLepton.Pt();
          gt->genLep2Eta   = dressedLepton.Eta();
          gt->genLep2Phi   = dressedLepton.Phi();
          gt->genLep2PdgId = thePdgId; 
        }
 
        if(v1.Pt() > 0 && DeltaR2(part.eta(),part.phi(),v1.Eta(),v1.Phi()) < 0.01) {
	  if     (part.testFlag(GenParticle::kIsTauDecayProduct) || part.testFlag(GenParticle::kIsPromptTauDecayProduct) || 
	          part.testFlag(GenParticle::kIsDirectTauDecayProduct) || part.testFlag(GenParticle::kIsDirectPromptTauDecayProduct) ||
		  (part.parent.isValid() && abs(part.parent->pdgid) == 15)) gt->looseGenLep1PdgId = 2;
	  else if(part.testFlag(GenParticle::kIsPrompt) || part.statusFlags == GenParticle::kIsPrompt) gt->looseGenLep1PdgId = 1;
	  if(part.pdgid != gt->looseLep1PdgId) gt->looseGenLep1PdgId = -1 * gt->looseGenLep1PdgId;
	}

        if(v2.Pt() > 0 && DeltaR2(part.eta(),part.phi(),v2.Eta(),v2.Phi()) < 0.01) {
	  if     (part.testFlag(GenParticle::kIsTauDecayProduct) || part.testFlag(GenParticle::kIsPromptTauDecayProduct) || 
	          part.testFlag(GenParticle::kIsDirectTauDecayProduct) || part.testFlag(GenParticle::kIsDirectPromptTauDecayProduct) ||
		  (part.parent.isValid() && abs(part.parent->pdgid) == 15)) gt->looseGenLep2PdgId = 2;
	  else if(part.testFlag(GenParticle::kIsPrompt) || part.statusFlags == GenParticle::kIsPrompt) gt->looseGenLep2PdgId = 1;
	  if(part.pdgid != gt->looseLep2PdgId) gt->looseGenLep2PdgId = -1 * gt->looseGenLep2PdgId;
	}

        if(v3.Pt() > 0 && DeltaR2(part.eta(),part.phi(),v3.Eta(),v3.Phi()) < 0.01) {
	  if     (part.testFlag(GenParticle::kIsTauDecayProduct) || part.testFlag(GenParticle::kIsPromptTauDecayProduct) || 
	          part.testFlag(GenParticle::kIsDirectTauDecayProduct) || part.testFlag(GenParticle::kIsDirectPromptTauDecayProduct) ||
		  (part.parent.isValid() && abs(part.parent->pdgid) == 15)) gt->looseGenLep3PdgId = 2;
	  else if(part.testFlag(GenParticle::kIsPrompt) || part.statusFlags == GenParticle::kIsPrompt) gt->looseGenLep3PdgId = 1;
	  if(part.pdgid != gt->looseLep3PdgId) gt->looseGenLep3PdgId = -1 * gt->looseGenLep3PdgId;
	}

        if(v4.Pt() > 0 && DeltaR2(part.eta(),part.phi(),v4.Eta(),v4.Phi()) < 0.01) {
	  if     (part.testFlag(GenParticle::kIsTauDecayProduct) || part.testFlag(GenParticle::kIsPromptTauDecayProduct) || 
	          part.testFlag(GenParticle::kIsDirectTauDecayProduct) || part.testFlag(GenParticle::kIsDirectPromptTauDecayProduct) ||
		  (part.parent.isValid() && abs(part.parent->pdgid) == 15)) gt->looseGenLep4PdgId = 2;
	  else if(part.testFlag(GenParticle::kIsPrompt) || part.statusFlags == GenParticle::kIsPrompt) gt->looseGenLep4PdgId = 1;
	  if(part.pdgid != gt->looseLep4PdgId) gt->looseGenLep4PdgId = -1 * gt->looseGenLep4PdgId;
	}
      }
      
      // Filling Z info at gen level
      if(gt->genLep1Pt > 25 && TMath::Abs(gt->genLep1Eta) < 2.4 && abs(gt->genLep1PdgId) != 15 &&
         gt->genLep2Pt > 25 && TMath::Abs(gt->genLep2Eta) < 2.4 && abs(gt->genLep2PdgId) != 15){
        TLorentzVector genlep1;
        genlep1.SetPtEtaPhiM(gt->genLep1Pt,gt->genLep1Eta,gt->genLep1Phi,0.0);
        TLorentzVector genlep2;
        genlep2.SetPtEtaPhiM(gt->genLep2Pt,gt->genLep2Eta,gt->genLep2Phi,0.0);
	TLorentzVector dilep = genlep1 + genlep2;
	if(TMath::Abs(dilep.M()-91.1876) < 15.0) {
          double ZGenPt  = TMath::Min(dilep.Pt(),1499.999);
          double ZGenRap = TMath::Abs(dilep.Rapidity());
	  double the_phi_star_eta = TMath::Min(phi_star_eta(genlep1,genlep2,gt->looseLep1PdgId),49.999);
	  if(the_phi_star_eta <= 1e-3) the_phi_star_eta = 1e-3;

	  if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtMM->Fill(ZGenPt,event.weight);
	  else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtEE->Fill(ZGenPt,event.weight);
	  if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPhiStarMM->Fill(the_phi_star_eta,event.weight);
	  else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPhiStarEE->Fill(the_phi_star_eta,event.weight);

	  if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtMM_PDF->Fill(ZGenPt,event.weight*gt->pdfUp);
	  else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtEE_PDF->Fill(ZGenPt,event.weight*gt->pdfUp);
	  if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPhiStarMM_PDF->Fill(the_phi_star_eta,event.weight*gt->pdfUp);
	  else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPhiStarEE_PDF->Fill(the_phi_star_eta,event.weight*gt->pdfUp);

          for(int i=0; i<6; i++){
	    if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtMM_QCDPart[i]->Fill(ZGenPt,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
	    else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtEE_QCDPart[i]->Fill(ZGenPt,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
	    if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPhiStarMM_QCDPart[i]->Fill(the_phi_star_eta,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
	    else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPhiStarEE_QCDPart[i]->Fill(the_phi_star_eta,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
          }

	  if(ZGenRap < 2.4) {
	    if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilRapMM->Fill(ZGenRap,event.weight);
	    else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilRapEE->Fill(ZGenRap,event.weight);

	    if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilRapMM_PDF->Fill(ZGenRap,event.weight*gt->pdfUp);
	    else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilRapEE_PDF->Fill(ZGenRap,event.weight*gt->pdfUp);

            for(int i=0; i<6; i++){
	      if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilRapMM_QCDPart[i]->Fill(ZGenRap,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
	      else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilRapEE_QCDPart[i]->Fill(ZGenRap,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
            }
	  }
	  if     (ZGenRap < 0.4) {
	    if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap0MM->Fill(ZGenPt,event.weight);
	    else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap0EE->Fill(ZGenPt,event.weight);

	    if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap0MM_PDF->Fill(ZGenPt,event.weight*gt->pdfUp);
	    else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap0EE_PDF->Fill(ZGenPt,event.weight*gt->pdfUp);

            for(int i=0; i<6; i++){
	      if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap0MM_QCDPart[i]->Fill(ZGenPt,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
	      else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap0EE_QCDPart[i]->Fill(ZGenPt,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
            }
	  }
	  else if(ZGenRap < 0.8) {
	    if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap1MM->Fill(ZGenPt,event.weight);
	    else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap1EE->Fill(ZGenPt,event.weight);

	    if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap1MM_PDF->Fill(ZGenPt,event.weight*gt->pdfUp);
	    else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap1EE_PDF->Fill(ZGenPt,event.weight*gt->pdfUp);

            for(int i=0; i<6; i++){
	      if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap1MM_QCDPart[i]->Fill(ZGenPt,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
	      else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap1EE_QCDPart[i]->Fill(ZGenPt,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
            }
	  }
	  else if(ZGenRap < 1.2) {
	    if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap2MM->Fill(ZGenPt,event.weight);
	    else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap2EE->Fill(ZGenPt,event.weight);

	    if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap2MM_PDF->Fill(ZGenPt,event.weight*gt->pdfUp);
	    else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap2EE_PDF->Fill(ZGenPt,event.weight*gt->pdfUp);

            for(int i=0; i<6; i++){
	      if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap2MM_QCDPart[i]->Fill(ZGenPt,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
	      else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap2EE_QCDPart[i]->Fill(ZGenPt,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
            }
	  }
	  else if(ZGenRap < 1.6) {
	    if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap3MM->Fill(ZGenPt,event.weight);
	    else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap3EE->Fill(ZGenPt,event.weight);

	    if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap3MM_PDF->Fill(ZGenPt,event.weight*gt->pdfUp);
	    else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap3EE_PDF->Fill(ZGenPt,event.weight*gt->pdfUp);

            for(int i=0; i<6; i++){
	      if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap3MM_QCDPart[i]->Fill(ZGenPt,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
	      else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap3EE_QCDPart[i]->Fill(ZGenPt,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
            }
	  }
	  else if(ZGenRap < 2.4) {
	    if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap4MM->Fill(ZGenPt,event.weight);
	    else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap4EE->Fill(ZGenPt,event.weight);

	    if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap4MM_PDF->Fill(ZGenPt,event.weight*gt->pdfUp);
	    else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap4EE_PDF->Fill(ZGenPt,event.weight*gt->pdfUp);

            for(int i=0; i<6; i++){
	      if     (TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 13) hDDilPtRap4MM_QCDPart[i]->Fill(ZGenPt,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
	      else if(TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 11) hDDilPtRap4EE_QCDPart[i]->Fill(ZGenPt,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
            }
	  }
	}
      }

      // Filling WWSS info at gen level
      if(gt->genLep1Pt > 20 && TMath::Abs(gt->genLep1Eta) < 2.5 && abs(gt->genLep1PdgId) != 15 && 
         gt->genLep2Pt > 20 && TMath::Abs(gt->genLep2Eta) < 2.5 && abs(gt->genLep2PdgId) != 15 &&
	 deltaEtaJJGen > 2.5 && mJJGen > 500){
        TLorentzVector genlep1;
        genlep1.SetPtEtaPhiM(gt->genLep1Pt,gt->genLep1Eta,gt->genLep1Phi,0.0);
        TLorentzVector genlep2;
        genlep2.SetPtEtaPhiM(gt->genLep2Pt,gt->genLep2Eta,gt->genLep2Phi,0.0);
	TLorentzVector dilep = genlep1 + genlep2;
        hDWWSSMLL->Fill(TMath::Min((double)dilep.M(),599.999),event.weight);
      }

      // Filling WW info at gen level
      if(gt->genLep1Pt > 20 && TMath::Abs(gt->genLep1Eta) < 2.5 && 
         gt->genLep2Pt > 20 && TMath::Abs(gt->genLep2Eta) < 2.5 &&
	 ((TMath::Abs(gt->genLep1PdgId) == 13 && TMath::Abs(gt->genLep2PdgId) == 11) ||
	  (TMath::Abs(gt->genLep1PdgId) == 11 && TMath::Abs(gt->genLep2PdgId) == 13))){
        TLorentzVector genlep1;
        genlep1.SetPtEtaPhiM(gt->genLep1Pt,gt->genLep1Eta,gt->genLep1Phi,0.0);
        TLorentzVector genlep2;
        genlep2.SetPtEtaPhiM(gt->genLep2Pt,gt->genLep2Eta,gt->genLep2Phi,0.0);
	TLorentzVector dilep = genlep1 + genlep2;
	double mll    = TMath::Min((double)dilep.M(),799.999);
	double ptl1   = TMath::Min((double)gt->genLep1Pt,149.999);
	double ptl2   = TMath::Min((double)gt->genLep2Pt,149.999);
	double dphill = TMath::Abs(genlep1.DeltaPhi(genlep2));
	if(mll > 20.0) {
          hDWWMLL    ->Fill(mll,event.weight);
	  hDWWMLL_PDF->Fill(mll,event.weight*gt->pdfUp);
          for(int i=0; i<6; i++){
	    hDWWMLL_QCDPart[i]->Fill(mll,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
          }
        }
	if(ptl1 > 25.0) {
          hDWWPTL1    ->Fill(ptl1,event.weight);
	  hDWWPTL1_PDF->Fill(ptl1,event.weight*gt->pdfUp);
          for(int i=0; i<6; i++){
	    hDWWPTL1_QCDPart[i]->Fill(ptl1,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
          }
        }
	if(ptl2 > 20.0) {
          hDWWPTL2    ->Fill(ptl2,event.weight);
	  hDWWPTL2_PDF->Fill(ptl2,event.weight*gt->pdfUp);
          for(int i=0; i<6; i++){
	    hDWWPTL2_QCDPart[i]->Fill(ptl2,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
          }
          hDWWDPHILL    ->Fill(dphill,event.weight);
	  hDWWDPHILL_PDF->Fill(dphill,event.weight*gt->pdfUp);
          for(int i=0; i<6; i++){
	    hDWWDPHILL_QCDPart[i]->Fill(dphill,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
          }
        }
	if(nGoodGenJets[0] == 0) {
          hDWWN0JET    ->Fill(0.0,event.weight);
	  hDWWN0JET_PDF->Fill(0.0,event.weight*gt->pdfUp);
          for(int i=0; i<6; i++){
	    hDWWN0JET_QCDPart[i]->Fill(0.0,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
          }
	}
	if(nGoodGenJets[1] == 0) {
          hDWWN0JET    ->Fill(1.0,event.weight);
	  hDWWN0JET_PDF->Fill(1.0,event.weight*gt->pdfUp);
          for(int i=0; i<6; i++){
	    hDWWN0JET_QCDPart[i]->Fill(1.0,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
          }
	}
	if(nGoodGenJets[2] == 0) {
          hDWWN0JET    ->Fill(2.0,event.weight);
	  hDWWN0JET_PDF->Fill(2.0,event.weight*gt->pdfUp);
          for(int i=0; i<6; i++){
	    hDWWN0JET_QCDPart[i]->Fill(2.0,event.weight*TMath::Abs(1+gt->scale[i])/maxQCDscale);
          }
	}
      }

      for (int iG : targetsN) {
        auto& part(event.genParticles.at(iG));
        TLorentzVector neutrino;
        neutrino.SetPtEtaPhiM(part.pt(),part.eta(),part.phi(),part.m());

        // check there is no further copy:
        bool isLastCopy=true;
        for (int kG : targetsN) {
          if (event.genParticles.at(kG).parent.isValid() && event.genParticles.at(kG).parent.get() == &part) {
            isLastCopy=false;
            break;
          }
        }
        if (!isLastCopy)
          continue;
	
	the_rhoP4 = the_rhoP4 + neutrino;
      }

      TLorentzVector theZBosons(0,0,0,0);
      TLorentzVector theWBosons(0,0,0,0);
      int nZBosons = 0; int nWBosons = 0;
      for (int iG : targetsV) {
        auto& part(event.genParticles.at(iG));
        TLorentzVector boson;
        boson.SetPtEtaPhiM(part.pt(),part.eta(),part.phi(),part.m());

        // check there is no further copy:
        bool isLastCopy=true;
        for (int kG : targetsV) {
          if (event.genParticles.at(kG).parent.isValid() && event.genParticles.at(kG).parent.get() == &part) {
            isLastCopy=false;
            break;
          }
        }
        if (!isLastCopy)
          continue;
	
	if(boson.Pt() < bosonPtMin) bosonPtMin = boson.Pt();
	
	if(abs(part.pdgid) == 23) {theZBosons = theZBosons + boson; nZBosons++;}
	if(abs(part.pdgid) == 24) {theWBosons = theWBosons + boson; nWBosons++;}
      }
      if(nZBosons+nWBosons == 0) bosonPtMin = 0;

      if(nZBosons >= 2) {
        double the_rho = 0.0; if(the_rhoP4.P() > 0) the_rho = the_rhoP4.Pt()/the_rhoP4.P();
        double theZZCorr[2] {1,1};
        theZZCorr[0] = weightEWKCorr(bosonPtMin,1);
        float GENmZZ = theZBosons.M();
        theZZCorr[1] = kfactor_qqZZ_qcd_M(GENmZZ);
        gt->sf_zz = theZZCorr[0]*theZZCorr[1];
        if(the_rho <= 0.3) gt->sf_zzUnc = (1.0+TMath::Abs((theZZCorr[0]-1)*(15.99/9.89-1)));
	else               gt->sf_zzUnc = (1.0+TMath::Abs((theZZCorr[0]-1)               ));
      } else {
        gt->sf_zz    = 1.0;
	gt->sf_zzUnc = 1.0;
      }

      if(nWBosons == 1 && nZBosons == 1) {
        TLorentzVector theWZBoson = theWBosons + theZBosons;
        gt->sf_wz = weightEWKWZCorr(theWZBoson.M());
      } else {
        gt->sf_wz = 1.0;
      }
      
      if(nZBosons == 1) {
        gt->sf_zh     = weightZHEWKCorr(h1Corrs[ZHEwkCorr],     theZBosons.Pt());
        gt->sf_zhUp   = weightZHEWKCorr(h1Corrs[ZHEwkCorrUp],   theZBosons.Pt());
        gt->sf_zhDown = weightZHEWKCorr(h1Corrs[ZHEwkCorrDown], theZBosons.Pt());
      }
      else {
        gt->sf_zh     = 1.0;
        gt->sf_zhUp   = 1.0;
        gt->sf_zhDown = 1.0;
      }

      // ttbar pT weight
      TLorentzVector vT,vTbar;
      float pt_t=0, pt_tbar=0;
      for (int iG : targetsTop) {
        auto& part(event.genParticles.at(iG));

        // check there is no further copy:
        bool isLastCopy=true;
        for (int kG : targetsTop) {
          if (event.genParticles.at(kG).parent.isValid() && event.genParticles.at(kG).parent.get() == &part) {
            isLastCopy=false;
            break;
          }
        }
        if (!isLastCopy)
          continue;

        if (part.pdgid>0) {
         pt_t = part.pt();
        } else {
         pt_tbar = part.pt();
        }
     }
      if (pt_t>0 && pt_tbar>0) {
        gt->sf_tt = TMath::Sqrt(TMath::Exp(0.0615-0.0005*TMath::Min((float)400.,pt_t)) *
                  		TMath::Exp(0.0615-0.0005*TMath::Min((float)400.,pt_tbar)));
      }

    } // end gen study

    if (!PassPreselection())
      continue;

    tr.TriggerEvent("presel");

    if (!isData) {
      // now get the jet btag SFs
      vector<btagcand> btagcands;
      vector<double> sf_cent, sf_bUp, sf_bDown, sf_mUp, sf_mDown;
      vector<double> sf_cent_alt, sf_bUp_alt, sf_bDown_alt, sf_mUp_alt, sf_mDown_alt;

      unsigned int nJ = cleaned30Jets.size();
      for (unsigned int iJ=0; iJ!=nJ; ++iJ) {
        panda::Jet *jet = cleaned30Jets.at(iJ);
        int flavor=0;
        float genpt=0;
        for (auto& gen : event.genParticles) {
          int apdgid = abs(gen.pdgid);
          if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
            continue;
          double dr2 = DeltaR2(jet->eta(),jet->phi(),gen.eta(),gen.phi());
          if (dr2<0.09) {
            genpt = gen.pt();
            if (apdgid==4 || apdgid==5) {
              flavor=apdgid;
              break;
            } else {
              flavor=0;
            }
          }
        } // finding the jet flavor
        float pt = jet->pt();
        float btagUncFactor = 1;
        float eta = jet->eta();
        double eff(1),sf(1),sfUp(1),sfDown(1);
        unsigned int binpt = btagpt.bin(pt);
        unsigned int bineta = btageta.bin(fabs(eta));
        if (flavor==5)
          eff = beff[bineta][binpt];
        else if (flavor==4)
          eff = ceff[bineta][binpt];
        else
          eff = lfeff[bineta][binpt];
        if      (jet==cleaned30Jets.at(0)) {
          gt->jet1Flav = flavor;
          gt->jet1GenPt = genpt;
        } 
	else if (jet==cleaned30Jets.at(1)) {
          gt->jet2Flav = flavor;
          gt->jet2GenPt = genpt;
        }
	else if (jet==cleaned30Jets.at(2)) {
          gt->jet3Flav = flavor;
          gt->jet3GenPt = genpt;
        }
	else if (jet==cleaned30Jets.at(3)) {
          gt->jet4Flav = flavor;
          gt->jet4GenPt = genpt;
        }
      }
      // btagging study
      unsigned int nJ20 = cleaned20Jets.size();
      for (unsigned int iJ=0; iJ!=nJ20; ++iJ) {
        panda::Jet *jet = cleaned20Jets.at(iJ);
        // Need to repeat the operation since these are different jets
        int flavor=0;
        float genpt=0;
        for (auto& gen : event.genParticles) {
          int apdgid = abs(gen.pdgid);
          if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
            continue;
          double dr2 = DeltaR2(jet->eta(),jet->phi(),gen.eta(),gen.phi());
          if (dr2<0.09) {
            genpt = gen.pt();
            if (apdgid==4 || apdgid==5) {
              flavor=apdgid;
              break;
            } else {
              flavor=0;
            }
          }
        } // finding the jet flavor
        float pt = jet->pt();
        float btagUncFactor = 1;
        float eta = jet->eta();
        double eff(1),sf(1),sfUp(1),sfDown(1);
        unsigned int binpt = btagpt.bin(pt);
        unsigned int bineta = btageta.bin(fabs(eta));
        if (flavor==5)
          eff = beff[bineta][binpt];
        else if (flavor==4)
          eff = ceff[bineta][binpt];
        else
          eff = lfeff[bineta][binpt];

        CalcBJetSFs(bJetL,flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
        btagcands.push_back(btagcand(iJ,flavor,eff,sf,sfUp,sfDown));
        sf_cent.push_back(sf);

        if (flavor>0) {
          sf_bUp.push_back(sfUp); sf_bDown.push_back(sfDown);
          sf_mUp.push_back(sf); sf_mDown.push_back(sf);
        } else {
          sf_bUp.push_back(sf); sf_bDown.push_back(sf);
          sf_mUp.push_back(sfUp); sf_mDown.push_back(sfDown);
        }

      } // loop over jets

      EvalBTagSF(btagcands,sf_cent,GeneralLeptonicTree::bCent,GeneralLeptonicTree::bJet);
      EvalBTagSF(btagcands,sf_bUp,GeneralLeptonicTree::bBUp,GeneralLeptonicTree::bJet);
      EvalBTagSF(btagcands,sf_bDown,GeneralLeptonicTree::bBDown,GeneralLeptonicTree::bJet);
      EvalBTagSF(btagcands,sf_mUp,GeneralLeptonicTree::bMUp,GeneralLeptonicTree::bJet);
      EvalBTagSF(btagcands,sf_mDown,GeneralLeptonicTree::bMDown,GeneralLeptonicTree::bJet);

    }

    tr.TriggerEvent("ak4 gen-matching");

    gt->Fill();

  } // entry loop

  if (DEBUG) { PDebug("PandaLeptonicAnalyzer::Run","Done with entry loop"); }

} // Run()
