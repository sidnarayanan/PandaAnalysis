#include "../interface/PandaAnalyzer.h"
#include "TMath.h"
#include <algorithm>
#include <vector>


using namespace panda;
using namespace std;
using namespace pa;

#define ADDOP(X) ops_all.back()->addSubOp<X>()


PandaAnalyzer::PandaAnalyzer(Analysis* a, int debug_/*=0*/) :
  Analyzer("PandaAnalyzer", a, debug_),
  cfgop(analysis, gt, DEBUG),
  wIDs(v_make_shared<TString>())
{
  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Calling constructor");

  Config& cfg = cfgop.cfg;
  Utils& utils = cfgop.utils;

  gblop = new GlobalOp(event, cfg, utils, gt);
  ops_all.emplace_back(gblop);

  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Adding AnalysisOps");

  // Define analyses - user should not touch below code. 
  preselop = new ContainerOp("pre-sel", event, cfg, utils, gt);
  ops_all.emplace_back(preselop);
  ADDOP(MapOp);
  if (analysis.unpackedGen)
    ADDOP(DeepGenOp<UnpackedGenParticle>);
  else
    ADDOP(DeepGenOp<GenParticle>);
  ADDOP(TriggerOp);
  ADDOP(SimpleLeptonOp);
  ADDOP(ComplicatedLeptonOp);
  ADDOP(SimplePhotonOp);
  ADDOP(ComplicatedPhotonOp);
  ADDOP(RecoilOp);
  ADDOP(FatJetReclusterOp<GeneralTree>);
  ADDOP(FatJetOp);
  ADDOP(JetOp);
  ADDOP(TauOp);
  ADDOP(VBFCatOp);

  postselop = new ContainerOp("post-sel", event, cfg, utils, gt);
  ops_all.emplace_back(postselop);
  ADDOP(GenPOp);
  ADDOP(JetFlavorOp); 
  ADDOP(HbbMiscOp);
  ADDOP(GenVHOp); 
  ADDOP(KinFitOp);
  ADDOP(InclusiveLeptonOp);
  ADDOP(SoftActivityOp);
  ADDOP(FatJetMatchingOp);
  ADDOP(BTagSFOp);
  ADDOP(BTagWeightOp);
  ADDOP(TriggerEffOp);
  ADDOP(GenStudyEWKOp);
  ADDOP(QCDUncOp);
  ADDOP(GenLepOp);
  ADDOP(GenJetNuOp);
  ADDOP(HFCountingOp);
  ADDOP(KFactorOp);
  ADDOP(ZvvHClassOp);

  for (auto& op : ops_all)
    op->print();

  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Reading inputs");
  // Read inputs
  getInput();

  event.setStatus(*tIn, {"!*"});
  event.setAddress(*tIn, cfgop.get_inputBranches());

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

  TTree* tW = static_cast<TTree*>(fIn->Get("weights"));
  if (tW && analysis.processType == kSignal) {
    if (tW->GetEntries()!=377 && tW->GetEntries()!=22) {
      logger.error("PandaAnalyzer::PandaAnalyzer",
          TString::Format("Reweighting failed because only found %u weights!",
                          unsigned(tW->GetEntries())));
      throw runtime_error("");
    }
    TString *id = new TString();
    tW->SetBranchAddress("id",&id);
    unsigned nW = tW->GetEntriesFast();
    for (unsigned iW=0; iW!=nW; ++iW) {
      tW->GetEntry(iW);
      wIDs->push_back(*id);
    }
  } else if (analysis.processType==kSignal) {
    logger.error("PandaAnalyzer::PandaAnalyzer","This is a signal file, but the weights are missing!");
    throw runtime_error("");
  }
  registry.publishConst("wIDs", wIDs);

  // Define outputs
  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Writing outputs");
  gt.is_monohiggs      = (analysis.monoh || analysis.hbb);
  gt.is_hbb            = analysis.hbb;
  gt.is_vh             = analysis.vqqhbb; 
  gt.is_vbf            = analysis.vbf;
  gt.is_fatjet         = (analysis.fatjet || analysis.deepGen);
  gt.is_leptonic       = analysis.complicatedLeptons;
  gt.is_photonic       = analysis.complicatedPhotons;
  gt.is_monotop        = !(analysis.monoh || analysis.hbb || analysis.vbf);
  gt.is_breg           = analysis.bjetRegTraining;
  gt.btagWeights       = analysis.btagWeights;
  gt.useCMVA           = analysis.useCMVA;
  for (auto& id : *wIDs)
    gt.signal_weights[id] = 1;

  makeOutput(); 

  fOut->WriteTObject(hDTotalMCWeight); delete hDTotalMCWeight; hDTotalMCWeight = nullptr;
  if (hDNPUWeight != nullptr) {
    fOut->WriteTObject(hDNPUWeight); delete hDNPUWeight; hDNPUWeight = nullptr;
  }

  event.rng.setSize(20);

  // read input data
  cfgop.readData(analysis.datapath);
  for (auto& op : ops_all)
    op->readData(analysis.datapath);

  if (DEBUG) logger.debug("PandaAnalyzer::PandaAnalyzer","Called constructor");
}


PandaAnalyzer::~PandaAnalyzer()
{
}

void PandaAnalyzer::AddGoodLumiRange(int run, int l0, int l1)
{
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) { // don't know about this run yet
    vector<LumiRange> newLumiList;
    newLumiList.emplace_back(l0,l1);
    goodLumis[run] = newLumiList;
  } else {
    run_->second.emplace_back(l0,l1);
  }
}


bool PandaAnalyzer::PassGoodLumis(int run, int lumi)
{
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) {
    // matched no run
    if (DEBUG)
      logger.debug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i",run));
    return false;
  }

  // found the run, now look for a lumi range
  for (auto &range : run_->second) {
    if (range.Contains(lumi)) {
      if (DEBUG)
        logger.debug("PandaAnalyzer::PassGoodLumis",TString::Format("Accepting run=%i, lumi=%i",run,lumi));
      return true;
    }
  }

  // matched no lumi range
  if (DEBUG)
    logger.debug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i, lumi=%i",run,lumi));
  return false;
}


bool PandaAnalyzer::PassPresel(Selection::Stage stage)
{
  if (selections.size() == 0)
    return true;

  bool pass = false;
  for (auto& s : selections) {
    if (s->anded())
      continue;
    if (DEBUG>1)
      logger.debug("PandaAnalyzer::PassPresel",s->get_name());
    if (s->accept(stage)) {
      pass = true;
      break;
    }
  }

  for (auto& s : selections) {
    if (s->anded()) {
      if (DEBUG>1)
        logger.debug("PandaAnalyzer::PassPresel",s->get_name());
      pass = pass && s->accept(stage);
    }
  }

  return pass;
}



void PandaAnalyzer::Reset()
{
  for (auto& op : ops_all)
    op->reset();

  Analyzer::Reset();
}



void PandaAnalyzer::Terminate()
{
  for (auto& op : ops_all)
    op->terminate();

  Analyzer::Terminate();
}


// run
void PandaAnalyzer::Run()
{

  // INITIALIZE --------------------------------------------------------------------------
  unsigned nZero, nEvents, iE=0;
  setupRun(nZero, nEvents); 

  for (auto& op : ops_all)
    op->initialize(registry);

  ProgressReporter pr("PandaAnalyzer::Run",&iE,&nEvents,100);
  TimeReporter& tr = cfgop.cfg.tr;
  tr.TriggerEvent("configuration"); 

  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    pr.Report();

    Reset();
    event.getEntry(*tIn,iE);
    tr.TriggerEvent(TString::Format("GetEntry %u",iE));

    gblop->execute();

    if (analysis.isData && !PassGoodLumis(gt.runNumber,gt.lumiNumber))
        continue;

    preselop->execute();

    if (!PassPresel(Selection::sReco))
      continue;

    postselop->execute();

    if (!PassPresel(Selection::sGen)) // only check gen presel here
      continue;

    gt.Fill();

    tr.TriggerEvent("fill");

  }

  pr.Done();

  tr.Summary();
  for (auto& s : selections)
    s->report();

  if (DEBUG) { logger.debug("PandaAnalyzer::Run","Done with entry loop"); }

} // Run()
