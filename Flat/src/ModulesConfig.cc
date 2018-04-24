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

  // create an RNG
  event.rng.setSize(20);

  ////////////////////////////////////////////////////////////////////// 
  // manipulate which branches to read
  event.setStatus(*t, {"!*"}); // turn everything off first

  TString jetname = (analysis->puppiJets) ? "puppi" : "chs";
  panda::utils::BranchList readlist({"runNumber", "lumiNumber", "eventNumber","weight"});
  readlist.setVerbosity(0);

  if (analysis->genOnly) {
    readlist += {"genParticles","genReweight","ak4GenJets","genMet","genParticlesU","electrons"};
  } else { 
    readlist += {"runNumber", "lumiNumber", "eventNumber", "rho", 
                 "isData", "npv", "npvTrue", "weight", "chsAK4Jets", 
                 "electrons", "muons", "taus", "photons", 
                 "pfMet", "caloMet", "puppiMet", "rawMet", "trkMet", 
                 "recoil","metFilters","trkMet"};

    if (analysis->ak8) {
      readlist += {jetname+"AK8Jets", "subjets", jetname+"AK8Subjets","Subjets"};
      if (analysis->hbb)
        readlist.push_back("ak8GenJets");
    } else if (analysis->fatjet) {
      readlist += {jetname+"CA15Jets", "subjets", jetname+"CA15Subjets","Subjets"};
      if (analysis->hbb)
        readlist.push_back("ca15GenJets");
    }
    if (analysis->recluster || analysis->bjetRegression || 
        analysis->deep || analysis->hbb || analysis->complicatedPhotons) {
      readlist.push_back("pfCandidates");
    }
    if (analysis->deepTracks || analysis->bjetRegression || analysis->hbb) {
      readlist += {"tracks","vertices"};
    }

    if (analysis->bjetRegression || analysis->deepSVs)
      readlist.push_back("secondaryVertices");

    if (isData || analysis->applyMCTriggers) {
      readlist.push_back("triggers");
    }

    if (!isData) {
      readlist += {"genParticles","genReweight","ak4GenJets","genMet"};
      if (analysis->hbb)
        readlist.push_back("partons"); 
    }
  }

  event.setAddress(*t, readlist); // pass the readlist so only the relevant branches are turned on
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
    registry.publish("wIDs", &wIDs);
  } else if (analysis->processType==kSignal) {
    PError("PandaAnalyzer::Init","This is a signal file, but the weights are missing!");
    return 2;
  }


  ////////////////////////////////////////////////////////////////////// 

  // manipulate the output tree
  if (isData) {
    std::vector<TString> droppable = {"mcWeight","scale","scaleUp",
                                      "trueGenBosonPt",
                                      "scaleDown","pdf.*","gen.*","sf_.*"};
    gt->RemoveBranches(droppable,{"sf_phoPurity"});
  }
  if (analysis->genOnly) {
    std::vector<TString> keepable = {"mcWeight","scale","scaleUp",
                                     "scaleDown","pdf.*","gen.*",
                                     "sf_tt.*","sf_qcdTT.*",
                                     "trueGenBosonPt","sf_qcd.*","sf_ewk.*",
                                     "nHF.*", "nB.*"};
    if (analysis->deepGen) {
      keepable.push_back("fj1ECFN.*");
      keepable.push_back("fjTau.*");
    }
    gt->RemoveBranches({".*"},keepable);
  }
  if (!analysis->fatjet && !analysis->ak8) {
    gt->RemoveBranches({"fj.*"});
  }
  if (analysis->complicatedLeptons) {
    gt->RemoveBranches({"genJet.*","puppiU.*","pfU.*","dphipfU.*","dphipuppi.*","jet.*"});
  }


  ////////////////////////////////////////////////////////////////////// 

  // give selectors access to the tree
  for (auto* s : selections)
    s->set_gt(gt);

  ////////////////////////////////////////////////////////////////////// 

  configMod = new ConfigMod(analysis, DEBUG);

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

  gt->is_monohiggs      = (analysis->monoh || analysis->hbb);
  gt->is_vbf            = analysis->vbf;
  gt->is_fatjet         = (analysis->fatjet || analysis->deepGen);
  gt->is_leptonic       = analysis->complicatedLeptons;
  gt->is_photonic       = analysis->complicatedPhotons;
  gt->is_monotop        = !(analysis->monoh || analysis->hbb || analysis->vbf);
  gt->btagWeights    = analysis->btagWeights;
  gt->useCMVA        = analysis->useCMVA;

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
    gt->signal_weights[id] = 1;

  // Build the input tree here 
  gt->WriteTree(tOut);

  if (DEBUG) PDebug("PandaAnalyzer::SetOutputFile","Created output in "+fOutPath);
}




void PandaAnalyzer::OpenCorrection(CorrectionType ct, TString fpath, TString hname, int dim) 
{
  fCorrs[ct] = TFile::Open(fpath);
  if (dim==1) 
    h1Corrs[ct] = new THCorr1(fCorrs[ct]->Get(hname));
  else if (dim==2)
    h2Corrs[ct] = new THCorr2(fCorrs[ct]->Get(hname));
  else 
    f1Corrs[ct] = new TF1Corr((TF1*)fCorrs[ct]->Get(hname));
}


void PandaAnalyzer::SetDataDir(const char *s) 
{
  TString dirPath(s);
  dirPath += "/";

  

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





