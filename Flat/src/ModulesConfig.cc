#include "../interface/PandaAnalyzer.h"
#include "TSystem.h"
#include <algorithm>
#include <vector>


// contains any configuration-related method

using namespace panda;
using namespace std;

int PandaAnalyzer::Init(TTree *t, TH1D *hweights, TTree *weightNames)
{
  if (DEBUG) PDebug("PandaAnalyzer::Init","Starting initialization");
  if (!t || !hweights) {
    PError("PandaAnalyzer::Init","Malformed input!");
    return 0;
  }
  tIn = t;

  ////////////////////////////////////////////////////////////////////// 
  // instantiate a configuration
  configMod = new ConfigMod(*analysis, gt, DEBUG);
  mods.push_back(configMod);

  ////////////////////////////////////////////////////////////////////// 
  // manipulate which branches to read
  event.setStatus(*t, {"!*"}); // turn everything off first
  event.setAddress(*t, configMod.get_inputBranches()); 
  if (DEBUG) PDebug("PandaAnalyzer::Init","Set addresses");

  ////////////////////////////////////////////////////////////////////// 
  // read MC weights
  hDTotalMCWeight = dynamic_cast<TH1D*>(hweights->Clone("hDTotalMCWeight"));
  hDTotalMCWeight->SetDirectory(0);

  if (weightNames && analysis->processType==kSignal) { // hack?
    if (weightNames->GetEntries()!=377 && weightNames->GetEntries()!=22) {
      PError("PandaAnalyzer::Init",
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
    registry.publishConst("wIDs", &wIDs);
  } else if (analysis->processType==kSignal) {
    PError("PandaAnalyzer::Init","This is a signal file, but the weights are missing!");
    return 2;
  }

  ////////////////////////////////////////////////////////////////////// 
  // give selectors access to the tree
  for (auto* s : selections)
    s->set_gt(&gt);

  ////////////////////////////////////////////////////////////////////// 
  // create an RNG
  event.rng.setSize(20);

  if (DEBUG) PDebug("PandaAnalyzer::Init","Finished configuration");

  return 0;
}


void PandaAnalyzer::SetOutputFile(TString fOutName) 
{
  fOutPath = fOutName;
  fOut = new TFile(fOutName,"RECREATE");
  fOut->cd();
  tOut = new TTree("events","events");

  fOut->WriteTObject(hDTotalMCWeight);    

  gt.is_monohiggs      = (analysis->monoh || analysis->hbb);
  gt.is_vbf            = analysis->vbf;
  gt.is_fatjet         = (analysis->fatjet || analysis->deepGen);
  gt.is_leptonic       = analysis->complicatedLeptons;
  gt.is_photonic       = analysis->complicatedPhotons;
  gt.is_monotop        = !(analysis->monoh || analysis->hbb || analysis->vbf);
  gt.btagWeights       = analysis->btagWeights;
  gt.useCMVA           = analysis->useCMVA;

  if (analysis->deep) {
    auxFilePath = fOutName.ReplaceAll(".root","_pf_%u.root");
    IncrementAuxFile();
  }

  if (analysis->deepGen) {
    auxFilePath = fOutName.ReplaceAll(".root","_gen_%u.root");
    IncrementGenAuxFile();
  }

  // fill the signal weights
  for (auto& id : wIDs) 
    gt.signal_weights[id] = 1;

  // Build the input tree here 
  gt.WriteTree(tOut);

  if (DEBUG) PDebug("PandaAnalyzer::SetOutputFile","Created output in "+fOutPath);
}

void PandaAnalyzer::SetDataDir(const char *s) 
{
  for (auto* mod : mods)
    mod->readData(s);
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

