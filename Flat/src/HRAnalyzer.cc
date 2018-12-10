#include "../interface/HRAnalyzer.h"
#include "TMath.h"
#include <algorithm>
#include <vector>


using namespace panda;
using namespace std;
using namespace pa;
using namespace fastjet;


HRAnalyzer::HRAnalyzer(Analysis* a, int debug_/*=0*/) :
  Analyzer("HRAnalyzer", a, debug_),
  cfg(analysis, DEBUG)
{
  if (DEBUG) logger.debug("HRAnalyzer::HRAnalyzer","Calling constructor");

  utils.fMSDcorr.reset(new TFile(analysis.datapath+"/puppiCorr.root"));
  utils.puppisd_corrGEN = (TF1*)utils.fMSDcorr->Get("puppiJECcorr_gen");
  utils.puppisd_corrRECO_cen = (TF1*)utils.fMSDcorr->Get("puppiJECcorr_reco_0eta1v3"); 
  utils.puppisd_corrRECO_for = (TF1*)utils.fMSDcorr->Get("puppiJECcorr_reco_1v3eta2v5"); 

  int activeAreaRepeats = 1;
  double ghostArea = 0.01;
  double ghostEtaMax = 7.0;
  utils.activeArea.reset(new GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea));
  utils.areaDef.reset(new AreaDefinition(active_area_explicit_ghosts,*(utils.activeArea)));
  double radius = 1.5;
  double sdZcut = 0.15;
  double sdBeta = 1.;
  if (analysis.ak8) {
    radius = 0.8;
    sdZcut = 0.1;
    sdBeta = 0.;
  }
  utils.softDrop.reset(new contrib::SoftDrop(sdBeta,sdZcut,radius));

  gen.reset(new HRGenPOp(event, cfg, utils, gt));
  fj.reset(new FatJetReclusterOp<HeavyResTree>(event, cfg, utils, gt));
  op.reset(new HRTagOp(event, cfg, utils, gt));

  gen->print();
  fj->print();
  op->print();

  if (DEBUG) logger.debug("HRAnalyzer::HRAnalyzer","Reading inputs");

  // Read inputs
  getInput();
  event.setStatus(*tIn, {"!*"});
  utils::BranchList bl;
  bl += {"runNumber", "lumiNumber", "eventNumber", "rho",
         "npv", "genParticles", "mcWeight", "recoil" , "pfCandidates"}; 
  if (analysis.ak8) {
    bl += {"puppiAK8Jets", "subjets", "puppiAK8Subjets","Subjets"};
  } else {
    bl += {"puppiCA15Jets", "subjets", "puppiCA15Subjets","Subjets"};
  }
  event.setAddress(*tIn, bl);

  TH1D* hDTotalMCWeight = static_cast<TH1D*>(static_cast<TH1D*>(fIn->Get("hSumW"))->Clone("hDTotalMCWeight"));
  hDTotalMCWeight->SetDirectory(0);
  TH1D* hDNPUWeight = nullptr; {
    TH1D* hbase = static_cast<TH1D*>(fIn->Get("hNPVTrue"));
    if (hbase == nullptr)
      hbase = static_cast<TH1D*>(fIn->Get("hNPVReco"));
    if (hbase != nullptr) {
      hDNPUWeight = static_cast<TH1D*>(hbase->Clone("hDNPUWeight"));
      hDNPUWeight->SetDirectory(0);
    }
  }

  // Define outputs

  if (DEBUG) logger.debug("HRAnalyzer::HRAnalyzer","Writing outputs");

  makeOutput(); 

  fOut->WriteTObject(hDTotalMCWeight); delete hDTotalMCWeight; hDTotalMCWeight = nullptr;
  if (hDNPUWeight != nullptr) {
    fOut->WriteTObject(hDNPUWeight); delete hDNPUWeight; hDNPUWeight = nullptr;
  }


  event.rng.setSize(20);

  // read input data
  gen->readData(analysis.datapath);
  fj->readData(analysis.datapath);
  op->readData(analysis.datapath);

  if (DEBUG) logger.debug("HRAnalyzer::HRAnalyzer","Called constructor");
}


HRAnalyzer::~HRAnalyzer()
{
}

void HRAnalyzer::Reset()
{
  gen->reset();
  fj->reset();
  op->reset();

  Analyzer::Reset();
}



void HRAnalyzer::Terminate()
{
  gen->terminate();
  fj->terminate();
  op->terminate();

  Analyzer::Terminate();
}


// run
void HRAnalyzer::Run()
{

  // INITIALIZE --------------------------------------------------------------------------

  unsigned nZero, nEvents, iE=0;
  setupRun(nZero, nEvents); 

  gen->initialize(registry);
  fj->initialize(registry);
  op->initialize(registry);

  ProgressReporter pr("HRAnalyzer::Run",&iE,&nEvents,100);
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

  if (DEBUG) { logger.debug("HRAnalyzer::Run","Done with entry loop"); }

} // Run()
