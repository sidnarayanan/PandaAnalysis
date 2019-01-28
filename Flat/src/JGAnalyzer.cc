#include "../interface/JGAnalyzer.h"
#include "TMath.h"
#include <algorithm>
#include <vector>


using namespace panda;
using namespace std;
using namespace pa;
using namespace fastjet;


JGAnalyzer::JGAnalyzer(Analysis* a, int debug_/*=0*/) :
  Analyzer("JGAnalyzer", a, debug_),
  cfg(analysis, DEBUG)
{
  if (DEBUG) logger.debug("JGAnalyzer::JGAnalyzer","Calling constructor");


  int activeAreaRepeats = 1;
  double ghostArea = 0.01;
  double ghostEtaMax = 7.0;
  utils.activeArea.reset(new GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea));
  utils.areaDef.reset(new AreaDefinition(active_area_explicit_ghosts,*(utils.activeArea)));

  gen.reset(new BaseGenPOp<JetGraphTree>(event, cfg, utils, gt));
  fj.reset(new FatJetReclusterOp<JetGraphTree>(event, cfg, utils, gt));
  op.reset(new JetGraphOp(event, cfg, utils, gt));

  gen->print();
  fj->print();
  op->print();

  if (DEBUG) logger.debug("JGAnalyzer::JGAnalyzer","Reading inputs");

  analysis.ak = true;
  analysis.ak8 = true;
  analysis.reclusterFJ = false;

  // Read inputs
  getInput();
  event.setStatus(*tIn, {"!*"});
  utils::BranchList bl;
  bl += {"runNumber", "lumiNumber", "eventNumber", "rho",
         "npv", "genParticles", "mcWeight", "recoil" , "pfCandidates", "puppiAK8Jets"}; 
  event.setAddress(*tIn, bl);

  TH1D* hDTotalMCWeight = static_cast<TH1D*>(static_cast<TH1D*>(fIn->Get("hSumW"))->Clone("hDTotalMCWeight"));
  hDTotalMCWeight->SetDirectory(0);

  // Define outputs

  if (DEBUG) logger.debug("JGAnalyzer::JGAnalyzer","Writing outputs");

  makeOutput(); 

  fOut->WriteTObject(hDTotalMCWeight); delete hDTotalMCWeight; hDTotalMCWeight = nullptr;


  event.rng.setSize(20);

  // read input data
  gen->readData(analysis.datapath);
  fj->readData(analysis.datapath);
  op->readData(analysis.datapath);

  if (DEBUG) logger.debug("JGAnalyzer::JGAnalyzer","Called constructor");
}


JGAnalyzer::~JGAnalyzer()
{
}

void JGAnalyzer::Reset()
{
  gen->reset();
  fj->reset();
  op->reset();

  Analyzer::Reset();
}



void JGAnalyzer::Terminate()
{
  gen->terminate();
  fj->terminate();
  op->terminate();

  Analyzer::Terminate();
}


// run
void JGAnalyzer::Run()
{

  // INITIALIZE --------------------------------------------------------------------------

  unsigned nZero, nEvents, iE=0;
  setupRun(nZero, nEvents); 

  gen->initialize(registry);
  fj->initialize(registry);
  op->initialize(registry);

  ProgressReporter pr("JGAnalyzer::Run",&iE,&nEvents,100);
  TimeReporter& tr = cfg.tr;

  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    tr.Start();
    pr.Report();

    Reset();
    event.getEntry(*tIn,iE);
    tr.TriggerEvent(TString::Format("GetEntry %u",iE));

    gen->execute();
    fj->execute();
    op->execute(); 
  }

  pr.Done(); 

  tr.Summary();

  if (DEBUG) { logger.debug("JGAnalyzer::Run","Done with entry loop"); }

} // Run()
