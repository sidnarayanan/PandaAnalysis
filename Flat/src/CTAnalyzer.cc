#include "../interface/CTAnalyzer.h"
#include "TMath.h"
#include <algorithm>
#include <vector>


using namespace std;
using namespace pa;


CTAnalyzer::CTAnalyzer(Analysis* a, int debug_/*=0*/) :
  Analyzer("CTAnalyzer", a, debug_)
{
  if (DEBUG) logger.debug("CTAnalyzer::CTAnalyzer","Calling constructor");

  op.reset(new CTOp(gt, event));

  if (DEBUG) logger.debug("CTAnalyzer::CTAnalyzer","Reading inputs");

  // Read inputs
  getInput("events");
  tIn->SetBranchAddress("is_signal_new", &event.is_signal_new);
  tIn->SetBranchAddress("E", &event.E);
  tIn->SetBranchAddress("PX", &event.PX);
  tIn->SetBranchAddress("PY", &event.PY);
  tIn->SetBranchAddress("PZ", &event.PZ);


  if (DEBUG) logger.debug("CTAnalyzer::CTAnalyzer","Writing outputs");
  makeOutput(); 

  // read input data
  if (DEBUG) logger.debug("CTAnalyzer::CTAnalyzer","Called constructor");
}


CTAnalyzer::~CTAnalyzer()
{
  if (DEBUG) logger.debug("CTAnalyzer::~CTAnalyzer","Calling destructor");

  if (DEBUG) logger.debug("CTAnalyzer::~CTAnalyzer","Called destructor");

}

void CTAnalyzer::Reset()
{
  op->reset();
  Analyzer::Reset();
}



void CTAnalyzer::Terminate()
{
  Analyzer::Terminate();
}


// run
void CTAnalyzer::Run()
{

  // INITIALIZE --------------------------------------------------------------------------

  unsigned nZero, nEvents, iE=0;
  setupRun(nZero, nEvents); 

  ProgressReporter pr("CTAnalyzer::Run",&iE,&nEvents,10);
  TimeReporter tr("CTAnalyzer", DEBUG);

  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    tr.Start();
    pr.Report();

    Reset();
    tIn->GetEntry(iE); 
    tr.TriggerEvent(TString::Format("GetEntry %u",iE));

    op->execute(); 
    tr.TriggerEvent("execute CT"); 
  }

  tr.Summary();

  if (DEBUG) { logger.debug("CTAnalyzer::Run","Done with entry loop"); }

} // Run()
