#include "../interface/TagAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>

using namespace panda;
using namespace std;

TagAnalyzer::TagAnalyzer(int debug_/*=0*/) {
  DEBUG = debug_;

  if (DEBUG) PDebug("TagAnalyzer::TagAnalyzer","Calling constructor");
  gt = new TagTree();
  if (DEBUG) PDebug("TagAnalyzer::TagAnalyzer","Built TagTree");
  ibetas = gt->get_ibetas();
  Ns = gt->get_Ns();
  orders = gt->get_orders();
  if (DEBUG) PDebug("TagAnalyzer::TagAnalyzer","Called constructor");
}


TagAnalyzer::~TagAnalyzer() {
  if (DEBUG) PDebug("TagAnalyzer::~TagAnalyzer","Calling destructor");
}


void TagAnalyzer::ResetBranches() {
  gt->Reset();
}


void TagAnalyzer::SetOutputFile(TString fOutName) {
  fOut = new TFile(fOutName,"RECREATE");
  tOut = new TTree("events","events");

  fOut->WriteTObject(hDTotalMCWeight);    

  // Build the input tree here 
  gt->WriteTree(tOut);

  if (DEBUG) PDebug("TagAnalyzer::SetOutputFile","Created output in "+fOutName);
}


int TagAnalyzer::Init(TTree *t, TH1D *hweights)
{
  if (DEBUG) PDebug("TagAnalyzer::Init","Starting initialization");
  if (!t || !hweights) {
    PError("TagAnalyzer::Init","Malformed input!");
    return 0;
  }
  tIn = t;

  event.setStatus(*t, {"!*"}); // turn everything off first

  panda::utils::BranchList readlist({"runNumber", "lumiNumber", "eventNumber", 
                                     "isData", "npv", "npvTrue", "weight", "chsAK4Jets", 
                                     "pfMet",  "puppiMet", 
                                     "recoil","metFilters",});
  readlist.setVerbosity(DEBUG);

  TString jetname = "puppi";
  readlist += {jetname+"CA15Jets", "subjets", jetname+"CA15Subjets","Subjets"};
  
  readlist.push_back("genParticles");
 
  event.setAddress(*t, readlist); // pass the readlist so only the relevant branches are turned on
 
  if (DEBUG) PDebug("TagAnalyzer::Init","Set addresses");

  hDTotalMCWeight = new TH1D("hDTotalMCWeight","hDTotalMCWeight",1,0,2);
  hDTotalMCWeight->SetBinContent(1,hweights->GetBinContent(1));

  if (DEBUG) PDebug("TagAnalyzer::Init","Finished configuration");

  return 0;
}


void TagAnalyzer::Terminate() {
  fOut->WriteTObject(tOut);
  fOut->Close();

  for (auto *f : fCorrs)
    if (f)
      f->Close();
  for (auto *h : h1Corrs)
    delete h;
  for (auto *h : h2Corrs)
    delete h;

  delete hDTotalMCWeight;
  if (DEBUG) PDebug("TagAnalyzer::Terminate","Finished with output");
}

void TagAnalyzer::OpenCorrection(CorrectionType ct, TString fpath, TString hname, int dim) {
  fCorrs[ct] = TFile::Open(fpath);
  if (dim==1) 
    h1Corrs[ct] = new THCorr1((TH1D*)fCorrs[ct]->Get(hname));
  else
    h2Corrs[ct] = new THCorr2((TH2D*)fCorrs[ct]->Get(hname));
}

double TagAnalyzer::GetCorr(CorrectionType ct, double x, double y) {
  if (h1Corrs[ct]!=0) {
    return h1Corrs[ct]->Eval(x); 
  } else if (h2Corrs[ct]!=0) {
    return h2Corrs[ct]->Eval(x,y);
  } else {
    PError("TagAnalyzer::GetCorr",
       TString::Format("No correction is defined for CorrectionType=%u",ct));
    return 1;
  }
}

void TagAnalyzer::SetDataDir(const char *s) {
  TString dirPath(s);
  dirPath += "/";

  if (DEBUG) PDebug("TagAnalyzer::SetDataDir","Starting loading of data");


}


bool TagAnalyzer::PassPreselection() {
  return gt->partonPt > 250;
}

bool TagAnalyzer::CheckParton(panda::GenParticle *p, double &size) {
  if (processType == kTop) {
    if (abs(p->pdgid) != 6)
      return false; 

    panda::GenParticle *q1 = NULL, *q2 = NULL, *b = NULL;
    bool notLastCopy = false;
    for (auto &child : event.genParticles) {
      unsigned apdgid = abs(child.pdgid);
      if (apdgid == 6 && child.parent.isValid() && 
          child.parent.get() ==  p) {
        notLastCopy = true;
        break;
      }
      if (apdgid > 5)
        continue;
      if (apdgid == 5) {
        if (child.parent.isValid() && child.parent.get() == p) {
          b = &child; 
        }
      } else {
        if (child.parent.isValid() && abs(child.parent->pdgid) == 24) {
          bool foundTopParent = false; // go up the history, see if this W came from the top in question
          auto *parent = child.parent.get(); // initialize to the W
          while (parent->parent.isValid()) {
            parent = parent->parent.get();
            if (parent == p) {
              foundTopParent = true; 
              break;
            }
          }
          if (foundTopParent) {
            if (q1)
              q2 = &child; 
            else
              q1 = &child;
          }
        }  
      }
    }
    if (notLastCopy)
      return false;
    if (!(b && q1 && q2))
      return false;
    size = DeltaR2(p->eta(),p->phi(),q1->eta(),q1->phi());
    size = std::max(size,DeltaR2(p->eta(),p->phi(),q2->eta(),q2->phi()));
    size = std::max(size,DeltaR2(p->eta(),p->phi(),b->eta(),b->phi()));
    return true;

  } else if (processType == kQCD) { // light-quark or gluon jet 
    unsigned apdgid = abs(p->pdgid);
    if (apdgid > 5 && apdgid != 22)
      return false;
    if ((p->statusFlags & (1<<panda::GenParticle::kIsHardProcess)) == 0)
      return false;
    size = 0.01;  // just some dummy
    return true;
  }

  return false;
}

panda::FatJet *TagAnalyzer::Match(panda::GenParticle *p, double distance) {
  for (auto &fj : event.puppiCA15Jets) {
    if (DeltaR2(fj.eta(),fj.phi(),p->eta(),p->phi()) < distance)
      return &fj;
  }
  return NULL;
}


// run
void TagAnalyzer::Run() {

  // INITIALIZE --------------------------------------------------------------------------

  unsigned int nEvents = tIn->GetEntries();
  unsigned int nZero = 0;
  if (lastEvent>=0 && lastEvent<(int)nEvents)
    nEvents = lastEvent;
  if (firstEvent>=0)
    nZero = firstEvent;

  if (!fOut || !tIn) {
    PError("TagAnalyzer::Run","NOT SETUP CORRECTLY");
    exit(1);
  }


  // panda::FatJetCollection* fatjets(0);
  // if (flags["fatjet"]) {
  //  if (flags["puppi"])
  //   fatjets = &event.puppiCA15Jets;
  //  else
  //   fatjets = &event.chsCA15Jets;
  // }

  // set up reporters
  unsigned int iE=0;
  ProgressReporter pr("TagAnalyzer::Run",&iE,&nEvents,10);
  TimeReporter tr("TagAnalyzer::Run",DEBUG);

  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    tr.Start();
    pr.Report();
    event.getEntry(*tIn,iE);


    tr.TriggerEvent(TString::Format("GetEntry %u",iE));
    if (DEBUG>2) {
      PDebug("TagAnalyzer::Run::Dump","");
      event.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("TagAnalyzer::Run::Dump","");
      event.photons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("TagAnalyzer::Run::Dump","");
      event.muons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("TagAnalyzer::Run::Dump","");
      event.electrons.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("TagAnalyzer::Run::Dump","");
      event.chsAK4Jets.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("TagAnalyzer::Run::Dump","");
      event.pfMet.print(std::cout, 2);
      std::cout << std::endl;
      PDebug("TagAnalyzer::Run::Dump","");
      event.metMuOnlyFix.print(std::cout, 2);
      std::cout << std::endl;
    }

    if (event.recoil.max<175) // no ECFs below this
      continue;

    tr.TriggerEvent("initialize");
    for (auto &gen : event.genParticles) {
      ResetBranches();
      
      // event info
      //gt->mcWeight = (event.weight>0) ? 1 : -1;
      gt->mcWeight = event.weight;
      gt->runNumber = event.runNumber;
      gt->lumiNumber = event.lumiNumber;
      gt->eventNumber = event.eventNumber;
      gt->npv = event.npv;
      gt->pu = event.npvTrue;

      // met
      gt->pfmet = event.pfMet.pt;
      gt->puppimet = event.puppiMet.pt;

      double size = -1;
      if (gen.pt() < 250)
        continue;
      if (!CheckParton(&gen,size))
        continue;
      gt->partonPt = gen.pt();
      gt->partonEta = gen.eta();
      gt->partonSize = size;
      gt->partonPdgId = gen.pdgid;

      auto *reco = Match(&gen,0.36);
      if (reco) {
        gt->partonIsReco = 1;
        auto &fj = *reco;

        gt->fj1Pt = fj.pt();
        gt->fj1Eta = fj.eta();
        gt->fj1Phi = fj.phi();
        gt->fj1M = fj.m();
        gt->fj1MSD = fj.mSD;
        gt->fj1RawPt = fj.rawPt;

        // now we do substructure
        gt->fj1Tau32 = clean(fj.tau3/fj.tau2);
        gt->fj1Tau32SD = clean(fj.tau3SD/fj.tau2SD);
        gt->fj1Tau21 = clean(fj.tau2/fj.tau1);
        gt->fj1Tau21SD = clean(fj.tau2SD/fj.tau1SD);

        for (auto ibeta : ibetas) {
          for (auto N : Ns) {
            for (auto order : orders) {
              TagTree::ECFParams p;
              p.order = order; p.N = N; p.ibeta = ibeta;
              if (gt->fj1IsClean || true)
                gt->fj1ECFNs[p] = fj.get_ecf(order,N,ibeta);
              else
                gt->fj1ECFNs[p] = fj.get_ecf(order,N,ibeta);
            }
          }
        } //loop over betas
        gt->fj1HTTMass = fj.htt_mass;
        gt->fj1HTTFRec = fj.htt_frec;

        std::vector<panda::MicroJet const*> subjets;
        for (unsigned iS(0); iS != fj.subjets.size(); ++iS)
         subjets.push_back(&fj.subjets.objAt(iS));

        auto csvsort([](panda::MicroJet const* j1, panda::MicroJet const* j2)->bool {
          return j1->csv > j2->csv;
         });

        std::sort(subjets.begin(),subjets.end(),csvsort);
        if (subjets.size()>0) {
          gt->fj1MaxCSV = subjets.at(0)->csv;
          gt->fj1MinCSV = subjets.back()->csv;
          if (subjets.size()>1) {
            gt->fj1SubMaxCSV = subjets.at(1)->csv;
          }
        }
        gt->fj1DoubleCSV = fj.double_sub;

      } else {
        gt->partonIsReco = 0;
      }

      gt->Fill();
      tr.TriggerEvent("Fill parton");
    }


  } // entry loop

  if (DEBUG) { PDebug("TagAnalyzer::Run","Done with entry loop"); }

} // Run()

