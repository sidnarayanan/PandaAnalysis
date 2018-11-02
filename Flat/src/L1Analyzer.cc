#include "../interface/L1Analyzer.h"
#include "TMath.h"
#include <algorithm>
#include <vector>


using namespace std;
using namespace pa;


L1Analyzer::L1Analyzer(Analysis* a, int debug_/*=0*/) :
  Analyzer("L1Analyzer", a, debug_)
{
  if (DEBUG) logger.debug("L1Analyzer::L1Analyzer","Calling constructor");

  op.reset(new L1Op(gt, event));

  if (DEBUG) logger.debug("L1Analyzer::L1Analyzer","Reading inputs");

  // Read inputs
  getInput("ntuple/tree");
  tIn->SetBranchAddress("metFilter", &event.metFilter);
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
  makeOutput(); 

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
  op->reset();
  Analyzer::Reset();
}



void L1Analyzer::Terminate()
{
  Analyzer::Terminate();
}


// run
void L1Analyzer::Run()
{

  // INITIALIZE --------------------------------------------------------------------------

  unsigned nZero, nEvents, iE=0;
  setupRun(nZero, nEvents); 

  ProgressReporter pr("L1Analyzer::Run",&iE,&nEvents,10);
  TimeReporter tr("L1Analyzer", DEBUG);

  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    tr.Start();
    pr.Report();

    Reset();
    tIn->GetEntry(iE); 
    tr.TriggerEvent(TString::Format("GetEntry %u",iE));

    op->execute(); 
    tr.TriggerEvent("execute L1"); 

    gt.Fill();
  }

  tr.Summary();

  if (DEBUG) { logger.debug("L1Analyzer::Run","Done with entry loop"); }

} // Run()
