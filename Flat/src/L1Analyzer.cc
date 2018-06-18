#include "../interface/L1Analyzer.h"
#include "TMath.h"
#include <algorithm>
#include <vector>


using namespace std;
using namespace pa;


L1Analyzer::L1Analyzer(Analysis* a, int debug_/*=0*/) :
  DEBUG(debug_), 
  analysis(*a)
{
  if (DEBUG) logger.debug("L1Analyzer::L1Analyzer","Calling constructor");

  mod.reset(new L1Mod(gt, event));

  if (DEBUG) logger.debug("L1Analyzer::L1Analyzer","Reading inputs");

  // Read inputs
  fIn.reset(TFile::Open(analysis.inpath));
  tIn = static_cast<TTree*>(fIn->Get("ntuple/tree"));
  tIn->SetBranchAddress("run", &event.run);
  tIn->SetBranchAddress("lumi", &event.lumi);
  tIn->SetBranchAddress("event", &event.event);
  tIn->SetBranchAddress("bunchCrossing", &event.bunchCrossing);
  tIn->SetBranchAddress("jet_p4", &event.jet_p4);
  tIn->SetBranchAddress("jet_neutralEmFrac", &event.jet_neutralEmFrac);
  tIn->SetBranchAddress("jet_neutralHadFrac", &event.jet_neutralHadFrac);
  tIn->SetBranchAddress("met_p4", &event.met_p4);
  tIn->SetBranchAddress("L1EG_bx", &event.L1EG_bx);
  tIn->SetBranchAddress("L1EG_p4", &event.L1EG_p4);
  tIn->SetBranchAddress("L1EG_iso", &event.L1EG_iso);
  tIn->SetBranchAddress("L1GtBx", &event.L1GtBx);


  // Define outputs
  if (DEBUG) logger.debug("L1Analyzer::L1Analyzer","Writing outputs");
  fOut.reset(TFile::Open(analysis.outpath, "RECREATE"));
  fOut->cd();
  tOut = new TTree("events", "events");
  registry.publish("fOut", fOut);
  gt.WriteTree(tOut);

  // read input data
  if (DEBUG) logger.debug("L1Analyzer::L1Analyzer","Called constructor");
}


L1Analyzer::~L1Analyzer()
{
  if (DEBUG) logger.debug("L1Analyzer::~L1Analyzer","Calling destructor");

  if (DEBUG) logger.debug("L1Analyzer::~L1Analyzer","Called destructor");

}

void L1Analyzer::Reset()
{
  gt.Reset();
  mod->reset();
  if (DEBUG) logger.debug("L1Analyzer::Reset","Reset");
}



void L1Analyzer::Terminate()
{
  fOut->WriteTObject(tOut);
  fOut->Close();
  fOut = 0; tOut = 0;

  if (DEBUG) logger.debug("L1Analyzer::Terminate","Finished with output");
}


// run
void L1Analyzer::Run()
{

  // INITIALIZE --------------------------------------------------------------------------

  unsigned int nEvents = tIn->GetEntries();
  unsigned int nZero = 0;
  if (lastEvent >= 0 && lastEvent < (int)nEvents)
    nEvents = lastEvent;
  if (firstEvent >= 0)
    nZero = firstEvent;

  if (!fOut || !tIn) {
    logger.error("L1Analyzer::Run","NOT SETUP CORRECTLY");
    exit(1);
  }
  unsigned int iE=0;

  fOut->cd(); // to be absolutely sure

  ProgressReporter pr("L1Analyzer::Run",&iE,&nEvents,10);
  TimeReporter tr("L1Analyzer", DEBUG);

  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    tr.Start();
    pr.Report();

    Reset();
    tIn->GetEntry(iE); 
    tr.TriggerEvent(TString::Format("GetEntry %u",iE));

    mod->execute(); 
    tr.TriggerEvent("execute L1"); 

    gt.Fill();
  }

  tr.Summary();

  if (DEBUG) { logger.debug("L1Analyzer::Run","Done with entry loop"); }

} // Run()
