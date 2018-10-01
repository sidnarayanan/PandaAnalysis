#ifndef Analyzer_h
#define Analyzer_h

#include "AnalyzerUtilities.h"
#include "GeneralTree.h"
#include "Selection.h"
#include "Module.h"
#include "PandaAnalysis/Flat/interface/Common.h"


namespace pa {

  template<typename TREE>
  class Analyzer {
  public :
    Analyzer(TString name_, Analysis* a, int debug_=0) : 
      name(name_),
      analysis(*a),
      DEBUG(debug_) { }

    virtual ~Analyzer() {
      if (DEBUG) logger.debug("PandaAnalyzer::~PandaAnalyzer","Calling destructor");
      fIn->Close();
      if (DEBUG) logger.debug("PandaAnalyzer::~PandaAnalyzer","Called destructor");
    }

    Analyzer(const Analyzer&) = delete;
    Analyzer& operator=(const Analyzer&) = delete; 

    virtual void Reset() {
      gt.Reset();
      if (DEBUG) logger.debug(name+"::Reset","Reset");
    }

    virtual void Run() = 0;

    virtual void Terminate() {
      fOut->WriteTObject(tOut);
      fOut->Close();
      fOut = 0; tOut = 0;

      if (DEBUG) logger.debug(name+"::Terminate","Finished with output");
    }

    int firstEvent{0}, lastEvent{-1};

  protected:
    void makeOutput() {
      fOut.reset(TFile::Open(analysis.outpath, "RECREATE"));
      fOut->cd();
      tOut = new TTree("events", "events");
      gt.WriteTree(tOut);
      registry.publish("fOut", fOut);
    }

    void getInput(const char* treeName="events") {
      fIn.reset(TFile::Open(analysis.inpath));      
      tIn = static_cast<TTree*>(fIn->Get(treeName));
    }

    void setupRun(unsigned &nZero, unsigned &nEvents) {
      nEvents = tIn->GetEntries();
      nZero = 0;
      if (lastEvent >= 0 && lastEvent < (int)nEvents)
        nEvents = lastEvent;
      if (firstEvent >= 0)
        nZero = firstEvent;

      if (!fOut || !tIn) {
        logger.error(name+"::Run","NOT SETUP CORRECTLY");
        exit(1);
      }
      fOut->cd(); // to be absolutely sure
    }

    //////////////////////////////////////////////////////////////////////////////////////
    TString name; 
    TREE gt;
    Analysis& analysis; //!< configure what to run
    int DEBUG; //!< debug verbosity level
    Registry registry;

    // IO for the analyzer
    // output file is owned by Analyzer, but we use it elsewhere
    std::shared_ptr<TFile> fOut{nullptr};
    TTree *tOut{nullptr};

    std::unique_ptr<TFile> fIn{nullptr};
    TTree *tIn{nullptr};  // input tree to read

  };
}

#endif
