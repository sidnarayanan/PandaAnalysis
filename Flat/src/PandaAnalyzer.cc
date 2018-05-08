#include "../interface/PandaAnalyzer.h"
#include "TMath.h"
#include <algorithm>
#include <vector>


using namespace panda;
using namespace std;
using namespace pa; 

#define ADDRECO(X) mods_reco.push_back(new X(event, cfg, utils, gt)); 
#define ADDGEN(X) mods_gen.push_back(new X(event, cfg, utils, gt)); 


PandaAnalyzer::PandaAnalyzer(Analysis* a, int debug_/*=0*/) : 
  DEBUG(debug_),
  analysis(*a), 
  cfgmod(analysis, gt, DEBUG)
{
  if (DEBUG) PDebug("PandaAnalyzer::PandaAnalyzer","Calling constructor");
 
  Config& cfg = cfgmod.cfg;
  Utils& utils = cfgmod.utils; 

  gblmod = new GlobalMod(event, cfg, utils, gt);
  mods_all.push_back(gblmod);

  // Define analyses 
  ADDRECO(GenPMod)
  if (analysis.unpackedGen)
    ADDRECO(DeepGenMod<UnpackedGenParticle>)
  else
    ADDRECO(DeepGenMod<GenParticle>)
  ADDRECO(SimpleLeptonMod)
  ADDRECO(ComplicatedLeptonMod)
  ADDRECO(SimplePhotonMod)
  ADDRECO(ComplicatedPhotonMod)
  ADDRECO(RecoilMod)
  ADDRECO(FatJetMod)
  ADDRECO(JetMod)
  ADDRECO(TauMod)

  ADDGEN(InclusiveLeptonMod)
  ADDGEN(SoftActivityMod)
  ADDGEN(FatJetMatchingMod)
  ADDGEN(BTagSFMod)
  ADDGEN(BTagWeightMod)
  ADDGEN(TriggerEffMod)
  ADDGEN(GenStudyEWKMod)
  ADDGEN(QCDUncMod)
  ADDGEN(GenLepMod)
  ADDGEN(GenJetNuMod)
  ADDGEN(HFCountingMod)
  ADDGEN(KFactorMod)

  mods_all.insert(mods_all.end(), mods_reco.begin(), mods_reco.end());
  mods_all.insert(mods_all.end(), mods_gen.begin(), mods_gen.end());
  

  // Read inputs 
  fIn = TFile::Open(analysis.inpath);
  tIn = static_cast<TTree*>(fIn->Get("events"));
  event.setStatus(*tIn, {"!*"});
  event.setAddress(*tIn, cfgmod.get_inputBranches());

  TH1D* hDTotalMCWeight = static_cast<TH1D*>(static_cast<TH1D*>(fIn->Get("hSumW"))->Clone("hDTotalMCWeight"));
  hDTotalMCWeight->SetDirectory(0);
  TH1D* hDNPUWeight = nullptr; {
    TH1D* hbase = static_cast<TH1D*>(fIn->Get("hNPVTrue"));
    if (hbase != nullptr) {
      hDNPUWeight = static_cast<TH1D*>(hbase->Clone("hDNPUWeight"));
      hDNPUWeight->SetDirectory(0);
    }
  }

  TTree* tW = static_cast<TTree*>(fIn->Get("weights"));
  if (tW && analysis.processType == kSignal) {
    if (tW->GetEntries()!=377 && tW->GetEntries()!=22) {
      PError("PandaAnalyzer::PandaAnalyzer",
          TString::Format("Reweighting failed because only found %u weights!",
                          unsigned(tW->GetEntries())));
      throw runtime_error("");
    }
    TString *id = new TString();
    tW->SetBranchAddress("id",&id);
    unsigned nW = tW->GetEntriesFast();
    for (unsigned iW=0; iW!=nW; ++iW) {
      tW->GetEntry(iW);
      wIDs.push_back(*id);
    }
  } else if (analysis.processType==kSignal) {
    PError("PandaAnalyzer::PandaAnalyzer","This is a signal file, but the weights are missing!");
    throw runtime_error("");
  }
  registry.publishConst("wIDs", &wIDs);

  // Define outputs
  
  gt.is_monohiggs      = (analysis.monoh || analysis.hbb);
  gt.is_vbf            = analysis.vbf;
  gt.is_fatjet         = (analysis.fatjet || analysis.deepGen);
  gt.is_leptonic       = analysis.complicatedLeptons;
  gt.is_photonic       = analysis.complicatedPhotons;
  gt.is_monotop        = !(analysis.monoh || analysis.hbb || analysis.vbf);
  gt.btagWeights       = analysis.btagWeights;
  gt.useCMVA           = analysis.useCMVA;
  for (auto& id : wIDs)
    gt.signal_weights[id] = 1; 

  fOut = TFile::Open(analysis.outpath, "RECREATE");
  fOut->cd();
  tOut = new TTree("events", "events");

  fOut->WriteTObject(hDTotalMCWeight); delete hDTotalMCWeight; hDTotalMCWeight = nullptr; 
  if (hDNPUWeight != nullptr) {
    fOut->WriteTObject(hDNPUWeight); delete hDNPUWeight; hDNPUWeight = nullptr; 
  }

  registry.publish("fOut", fOut); 

  gt.WriteTree(tOut);

  event.rng.setSize(20);


  // read input data 
  cfgmod.readData(analysis.datapath); 
  for (auto* mod : mods_all)
    mod->readData(analysis.datapath);

  if (DEBUG) PDebug("PandaAnalyzer::PandaAnalyzer","Called constructor");
}


PandaAnalyzer::~PandaAnalyzer() 
{
  if (DEBUG) PDebug("PandaAnalyzer::~PandaAnalyzer","Calling destructor");
  
  for (auto* mod : mods_all)
    delete mod; 

  for (auto* s : selections)
    delete s; 

  fIn->Close();
  if (DEBUG) PDebug("PandaAnalyzer::~PandaAnalyzer","Called destructor");

}

void PandaAnalyzer::AddGoodLumiRange(int run, int l0, int l1) 
{
  auto run_ = goodLumis.find(run);
  if (run_==goodLumis.end()) { // don't know about this run yet
    std::vector<LumiRange> newLumiList;
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
      PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i",run));
    return false;
  }

  // found the run, now look for a lumi range
  for (auto &range : run_->second) {
    if (range.Contains(lumi)) {
      if (DEBUG) 
        PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Accepting run=%i, lumi=%i",run,lumi));
      return true;
    }
  }

  // matched no lumi range
  if (DEBUG) 
    PDebug("PandaAnalyzer::PassGoodLumis",TString::Format("Failing run=%i, lumi=%i",run,lumi));
  return false;
}


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



void PandaAnalyzer::Reset() 
{
  gt.Reset();
  
  for (auto* mod : mods_all)
    mod->reset(); 
  if (DEBUG) PDebug("PandaAnalyzer::Reset","Reset");
}



void PandaAnalyzer::Terminate() 
{
  fOut->WriteTObject(tOut);
  fOut->Close();
  fOut = 0; tOut = 0;

  for (auto* mod : mods_all)
    mod->terminate(); 

  if (DEBUG) PDebug("PandaAnalyzer::Terminate","Finished with output");
}


// run
void PandaAnalyzer::Run() 
{

  // INITIALIZE --------------------------------------------------------------------------

  unsigned int nEvents = tIn->GetEntries();
  unsigned int nZero = 0;
  if (lastEvent >= 0 && lastEvent < (int)nEvents)
    nEvents = lastEvent;
  if (firstEvent >= 0)
    nZero = firstEvent;

  if (!fOut || !tIn) {
    PError("PandaAnalyzer::Run","NOT SETUP CORRECTLY");
    exit(1);
  }
  unsigned int iE=0;

  fOut->cd(); // to be absolutely sure

  for (auto* mod : mods_all)
    mod->initialize(registry); 

  ProgressReporter pr("PandaAnalyzer::Run",&iE,&nEvents,10);
  TimeReporter& tr = cfgmod.cfg.tr; 

  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    tr.Start();
    pr.Report();

    Reset();
    event.getEntry(*tIn,iE);
    tr.TriggerEvent(TString::Format("GetEntry %u",iE));

    gblmod->execute();

    if (analysis.isData && !PassGoodLumis(gt.runNumber,gt.lumiNumber))
        continue;

    for (auto* mod : mods_reco)
      mod->execute();

    if (!PassPresel(Selection::sRecoil))
      continue;

    for (auto* mod : mods_gen)
      mod->execute();

    if (!PassPresel(Selection::sGen)) // only check gen presel here
      continue;

//     if (analysis.deep || analysis.hbb)
//       FatjetPartons();
//     if (analysis.deep) {
//       FillPFTree();
//       tAux->Fill();
//       if (tAux->GetEntriesFast() == 2500)
//         IncrementAuxFile();
//       tr->TriggerEvent("aux fill");
//     }

    gt.Fill();

    tr.TriggerEvent("fill");

  } 

  tr.Summary();
  for (auto* s : selections) 
    s->report(); 

  if (DEBUG) { PDebug("PandaAnalyzer::Run","Done with entry loop"); }

} // Run()

