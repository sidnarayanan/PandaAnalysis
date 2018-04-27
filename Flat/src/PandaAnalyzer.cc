#include "../interface/PandaAnalyzer.h"
#include "TMath.h"
#include <algorithm>
#include <vector>


// only contain the main methods

using namespace panda;
using namespace std;

PandaAnalyzer::PandaAnalyzer(int debug_/*=0*/) 
{
  DEBUG = debug_;

  if (DEBUG) PDebug("PandaAnalyzer::PandaAnalyzer","Calling constructor");
  if (DEBUG) PDebug("PandaAnalyzer::PandaAnalyzer","Called constructor");
}


PandaAnalyzer::~PandaAnalyzer() 
{
  if (DEBUG) PDebug("PandaAnalyzer::~PandaAnalyzer","Calling destructor");
}


void PandaAnalyzer::ResetBranches() 
{
  genObjects.clear();
  matchPhos.clear();
  matchEles.clear();
  matchLeps.clear();
  looseLeps.clear();
  tightLeps.clear();
  loosePhos.clear();
  genJetsNu.clear();
  validGenP.clear();
  for (int i = 0; i != jes2i(shiftjes::N); ++i) 
    jesShifts[i].clear();
  fj1 = 0;
  gt->Reset();
  if (DEBUG) PDebug("PandaAnalyzer::ResetBranches","Reset");
}



void PandaAnalyzer::Terminate() 
{
  fOut->WriteTObject(tOut);
  fOut->Close();
  fOut = 0; tOut = 0;

  if (analysis->deep)
    IncrementAuxFile(true);
  if (analysis->deepGen)
    IncrementGenAuxFile(true);

  if (DEBUG) PDebug("PandaAnalyzer::Terminate","Finished with output");
}


// run
void PandaAnalyzer::Run() 
{

  fOut->cd(); // to be absolutely sure

  // INITIALIZE --------------------------------------------------------------------------

  unsigned int nEvents = tIn->GetEntries();
  unsigned int nZero = 0;
  if (lastEvent>=0 && lastEvent<(int)nEvents)
    nEvents = lastEvent;
  if (firstEvent>=0)
    nZero = firstEvent;

  if (!fOut || !tIn) {
    PError("PandaAnalyzer::Run","NOT SETUP CORRECTLY");
    exit(1);
  }

  // get bounds
  genBosonPtMin=150, genBosonPtMax=1000;
  if (!isData && h1Corrs[cZNLO]) {
    genBosonPtMin = h1Corrs[cZNLO]->GetHist()->GetBinCenter(1);
    genBosonPtMax = h1Corrs[cZNLO]->GetHist()->GetBinCenter(h1Corrs[cZNLO]->GetHist()->GetNbinsX());
  }

  if (analysis->ak8) {
    if (analysis->puppiJets)
      fatjets = &event.puppiAK8Jets;
    else
      fatjets = &event.chsAK8Jets;
  } else if (analysis->fatjet) {
    if (analysis->puppiJets)
      fatjets = &event.puppiCA15Jets;
    else
      fatjets = &event.chsCA15Jets;
  }

  ak4jets = &event.chsAK4Jets;

  // these are bins of b-tagging eff in pT and eta, derived in 8024 TT MC
  // TODO: don't hardcode these 
  std::vector<double> vbtagpt {20.0,50.0,80.0,120.0,200.0,300.0,400.0,500.0,700.0,1000.0};
  std::vector<double> vbtageta {0.0,0.5,1.5,2.5};
  lfeff  = {{0.081,0.065,0.060,0.063,0.072,0.085,0.104,0.127,0.162},
            {0.116,0.097,0.092,0.099,0.112,0.138,0.166,0.185,0.222},
            {0.173,0.145,0.149,0.175,0.195,0.225,0.229,0.233,0.250}};
  ceff = {{0.377,0.389,0.391,0.390,0.391,0.375,0.372,0.392,0.435},
          {0.398,0.407,0.416,0.424,0.424,0.428,0.448,0.466,0.500},
          {0.375,0.389,0.400,0.425,0.437,0.459,0.481,0.534,0.488}};
  beff = {{0.791,0.815,0.825,0.835,0.821,0.799,0.784,0.767,0.760},
          {0.794,0.816,0.829,0.836,0.823,0.804,0.798,0.792,0.789},
          {0.739,0.767,0.780,0.789,0.776,0.771,0.779,0.787,0.806}};
  btagpt = Binner(vbtagpt);
  btageta = Binner(vbtageta);


  fOut->cd(); // to be absolutely sure

  // set up reporters
  unsigned int iE=0;
  ProgressReporter pr("PandaAnalyzer::Run",&iE,&nEvents,10);
  tr = new TimeReporter("PandaAnalyzer::Run",DEBUG+1);


  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    tr->Start();
    pr.Report();
    ResetBranches();
    event.getEntry(*tIn,iE);


    tr->TriggerEvent(TString::Format("GetEntry %u",iE));
    if (isData && !PassGoodLumis(gt.runNumber,gt.lumiNumber))
        continue;

    // global mod goes here

    if (!PassPresel(Selection::sRecoil))
      continue;


    tr->TriggerEvent("initialize");


    if (!isData) {
      // GenPMod

      // DeepGenMod<>
    }


    if (!analysis->genOnly) {
      // SimpleLepMod, ComplicatedLepMod
      // InclusiveLepMod
      //
      // SimplePhoMod, ComplicatedPhoMod
      //
      // RecoilMod
      //
      // FatJetMod 
      //   -> FatJetReclusterMod
      //
      // 

      // interesting jets
      JetBasics();

      Taus();

      if (!PassPresel(Selection::sReco)) // only check reco presel here
        continue;

      if (analysis->hbb) {
        JetHbbSoftActivity();
        GetMETSignificance();
      }
    }

    if (!isData) {
      if (!analysis->genOnly) {
        if (analysis->fatjet)
          FatjetMatching();

        if (analysis->btagSFs)
          JetBtagSFs();
        if (analysis->btagWeights)
          JetCMVAWeights();
        
        TriggerEffs();

        if (analysis->complicatedLeptons ||
            analysis->complicatedPhotons)
          GenStudyEWK();
        else
          LeptonSFs();

        PhotonSFs();
      }

      QCDUncs();
      if (analysis->hbb)
        LHEInfo();
      SignalReweights();

      if (analysis->vbf)
        SaveGenLeptons();

      SignalInfo();

      if (analysis->reclusterGen && analysis->hbb) {
        GenJetsNu();
        MatchGenJets(genJetsNu);
      }

      if (analysis->hfCounting)
        HeavyFlavorCounting();

       TopPTReweight();
       VJetsReweight();
    }

    
    if (!PassPresel(Selection::sGen)) // only check gen presel here
      continue;

    if (analysis->deep || analysis->hbb)
      FatjetPartons();
    if (analysis->deep) {
      FillPFTree();
      tAux->Fill();
      if (tAux->GetEntriesFast() == 2500)
        IncrementAuxFile();
      tr->TriggerEvent("aux fill");
    }

    gt->Fill();

    tr->TriggerEvent("fill");

  } // entry loop

  tr->Summary();
  for (auto* s : selections) 
    s->report(); 

  if (DEBUG) { PDebug("PandaAnalyzer::Run","Done with entry loop"); }

} // Run()

