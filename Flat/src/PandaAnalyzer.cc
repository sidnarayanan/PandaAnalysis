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

  fOut->cd(); // to be absolutely sure

  // set up reporters
  unsigned int iE=0;
  ProgressReporter pr("PandaAnalyzer::Run",&iE,&nEvents,10);


  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    cfg.tr->Start();
    pr.Report();
    ResetBranches();
    event.getEntry(*tIn,iE);


    cfg.tr->TriggerEvent(TString::Format("GetEntry %u",iE));
    if (analysis->isData && !PassGoodLumis(gt.runNumber,gt.lumiNumber))
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
      // JetMod
      //  -> JetFlavorMod
      //  -> IsoJetMod
      //  -> BJetRegMod
      //  -> VBFSystemMod
      //  -> HbbSystemMod 
      //
      // TauMod 

      if (!PassPresel(Selection::sReco)) // only check reco presel here
        continue;

      // SoftActivityMod
      // FatJetMatchingMod 
      // BTagSFMod
      // BTagWeightMod
      // TriggerEffMod
      // GenStudyEWKMod
      // QCDUncMod
      // GenLepMod
      // GenJetNuMod
      // HFCountingMod 
      // KFactorMod
    }

    if (!PassPresel(Selection::sGen)) // only check gen presel here
      continue;

//     if (analysis->deep || analysis->hbb)
//       FatjetPartons();
//     if (analysis->deep) {
//       FillPFTree();
//       tAux->Fill();
//       if (tAux->GetEntriesFast() == 2500)
//         IncrementAuxFile();
//       tr->TriggerEvent("aux fill");
//     }

    gt->Fill();

    tr->TriggerEvent("fill");

  } // entry loop

  tr->Summary();
  for (auto* s : selections) 
    s->report(); 

  if (DEBUG) { PDebug("PandaAnalyzer::Run","Done with entry loop"); }

} // Run()

