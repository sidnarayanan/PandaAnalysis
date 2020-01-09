#include "../interface/GrappleAnalyzer.h"
#include "TMath.h"
#include <algorithm>
#include <vector>


using namespace panda;
using namespace std;
using namespace pa;
using namespace fastjet;


GrappleAnalyzer::GrappleAnalyzer(Analysis* a, int debug_/*=0*/) :
  Analyzer("GrappleAnalyzer", a, debug_),
  cfg(analysis, DEBUG)
{
  if (DEBUG) logger.debug("GrappleAnalyzer::GrappleAnalyzer","Calling constructor");

  op.reset(new GrappleOp(event, cfg, utils, gt));

  op->print();

  cfg.NGENPROPS = 7;
  cfg.NMAXPF = 2000;
  cfg.auxFilePath = analysis.outpath;
  cfg.auxFilePath.ReplaceAll(".root", "_aux%i.root");

  if (DEBUG) logger.debug("GrappleAnalyzer::GrappleAnalyzer","Reading inputs");

  // Read inputs
  getInput();
  event.setStatus(*tIn, {"!*"});
  utils::BranchList bl;
  bl += {"runNumber", "lumiNumber", "eventNumber", "rho",
         "npv", "genParticles", "mcWeight", "pfCandidates", "vertices"}; 
  event.setAddress(*tIn, bl);

  TH1D* hDTotalMCWeight = static_cast<TH1D*>(static_cast<TH1D*>(fIn->Get("hSumW"))->Clone("hDTotalMCWeight"));
  hDTotalMCWeight->SetDirectory(0);

  // Define outputs

  if (DEBUG) logger.debug("GrappleAnalyzer::GrappleAnalyzer","Writing outputs");

  makeOutput(); 

  fOut->WriteTObject(hDTotalMCWeight); delete hDTotalMCWeight; hDTotalMCWeight = nullptr;

  event.rng.setSize(20);

  // read input data
  op->readData(analysis.datapath);

  if (DEBUG) logger.debug("GrappleAnalyzer::GrappleAnalyzer","Called constructor");
}


GrappleAnalyzer::~GrappleAnalyzer()
{
}

void GrappleAnalyzer::Reset()
{
  op->reset();

  Analyzer::Reset();
}



void GrappleAnalyzer::Terminate()
{
  op->terminate();

  Analyzer::Terminate();
}


// run
void GrappleAnalyzer::Run()
{

  // INITIALIZE --------------------------------------------------------------------------

  unsigned nZero, nEvents, iE=0;
  setupRun(nZero, nEvents); 

  op->initialize(registry);

  ProgressReporter pr("GrappleAnalyzer::Run",&iE,&nEvents,100);
  TimeReporter& tr = cfg.tr;

  // EVENTLOOP --------------------------------------------------------------------------
  for (iE=nZero; iE!=nEvents; ++iE) {
    tr.Start();
    pr.Report();

    Reset();
    event.getEntry(*tIn,iE);
    tr.TriggerEvent(TString::Format("GetEntry %u",iE));

    op->execute(); 
  }

  pr.Done(); 

  tr.Summary();

  if (DEBUG) { logger.debug("GrappleAnalyzer::Run","Done with entry loop"); }

} // Run()


void GrappleOp::incrementAux(bool close)
{
  if (fAux) {
    fAux->WriteTObject(tAux, "inputs", "Overwrite");
    TString path = TString::Format(cfg.auxFilePath.Data(),auxCounter);
    fAux->Close();
  }
  if (close)
    return;

  TString path = TString::Format(cfg.auxFilePath.Data(),auxCounter++);
  fAux = new TFile(path.Data(), "RECREATE");
  tAux = new TTree("inputs", "inputs");

  genJetInfo.particles.resize(cfg.NMAXPF);
  for (int i = 0; i != cfg.NMAXPF; ++i) {
    genJetInfo.particles[i].resize(cfg.NGENPROPS);
  }

  tAux->Branch("eventNumber",&(gt.eventNumber),"eventNumber/l");
  tAux->Branch("npv", &(gt.npv), "npv/I");
  tAux->Branch("kinematics",&(genJetInfo.particles));

  fOut->cd();
}


void GrappleOp::do_execute()
{

  gt.eventNumber = event.eventNumber;
  gt.npv = event.npv;

  std::map<const RecoVertex*, int> vertexToID;
  for (auto& vtx : event.vertices) {
    vertexToID[&vtx] = vertexToID.size();
  }

  auto sortedPFs = vector<const PFCand*>();
  for (auto& cand : event.pfCandidates)
    sortedPFs.push_back(&cand);
  sort(sortedPFs.begin(), sortedPFs.end(),
       [](const PFCand* x, const PFCand* y) -> bool {
          return x->pt() > y->pt();
        }); 

  int i_particle = -1;
  for (auto* candPtr : sortedPFs) {
    i_particle++;
    if (i_particle == cfg.NMAXPF)
      break;

    auto& cand = *candPtr;
    float pt = cand.pt();
    float eta = cand.eta();
    float phi = cand.phi();
    float e = cand.e();
    float puppi = cand.puppiW();

    // if particle is charged:
    //   vertexID = real vertex ID (0 if hard)
    // else:
    //   if particle is matched to hard scattering:
    //     vertexID = 0
    //   else:
    //     vertexID = -1
    int vertexID = -1;
    if (cand.q() != 0 && cand.vertex.isValid()) {
        vertexID = vertexToID[cand.vertex.get()];
    } else {
      for (auto& p : event.genParticles) {
        if (!p.finalState)
          continue;
        if (DeltaR2(eta, phi, p.eta(), p.phi()) < 0.0064) {
          vertexID = 0;
          break;
        }
      }   
    } 
    auto& entry = genJetInfo.particles[i_particle];
    entry[0] = pt;
    entry[1] = eta;
    entry[2] = phi;
    entry[3] = e;
    entry[4] = vertexID;
    entry[5] = cand.q();
    entry[6] = puppi;
  }

  tAux->Fill();
  if (tAux->GetEntries() == 10000)
    incrementAux(false);
}
