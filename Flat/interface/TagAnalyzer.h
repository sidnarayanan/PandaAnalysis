#ifndef TagAnalyzer_h
#define TagAnalyzer_h

// STL
#include "vector"
#include "map"
#include <string>
#include <cmath>

// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLorentzVector.h>

#include "AnalyzerUtilities.h"
#include "TagTree.h"


/////////////////////////////////////////////////////////////////////////////
// some misc definitions


/////////////////////////////////////////////////////////////////////////////
// TagAnalyzer definition
namespace pa {
    class TagAnalyzer {
    public :

        enum ProcessType { 
            kNone,
            kQCD,
            kTop,
        };


        //////////////////////////////////////////////////////////////////////////////////////

        TagAnalyzer(int debug_=0);
        ~TagAnalyzer();
        int Init(TTree *tree, TH1D *hweights);
        void SetOutputFile(TString fOutName);
        void ResetBranches();
        void Run();
        void Terminate();
        void SetDataDir(const char *s);

        // public configuration
        void SetFlag(TString flag, bool b=true) { flags[flag]=b; }
        int firstEvent=-1;
        int lastEvent=-1;                                                    // max events to process; -1=>all
        ProcessType processType=kNone;                         // determine what to do the jet matching to

    private:
        enum CorrectionType { //!< enum listing relevant corrections applied to MC
            cN=0
        };


        bool PassPreselection();
        void OpenCorrection(CorrectionType,TString,TString,int);
        double GetCorr(CorrectionType ct,double x, double y=0);
        bool CheckParton(panda::GenParticle*,double&);
        panda::FatJet *Match(panda::GenParticle*,double);

        int DEBUG = 0; //!< debug verbosity level
        std::map<TString,bool> flags;


        // files and histograms containing weights
        std::vector<TFile*> fCorrs = std::vector<TFile*>(cN,0); //!< files containing corrections
        std::vector<THCorr1*> h1Corrs = std::vector<THCorr1*>(cN,0); //!< histograms for binned corrections
        std::vector<THCorr2*> h2Corrs = std::vector<THCorr2*>(cN,0); //!< histograms for binned corrections

        // IO for the analyzer
        TFile *fOut;     // output file is owned by TagAnalyzer
        TTree *tOut;
        TagTree *gt; // essentially a wrapper around tOut
        TH1D *hDTotalMCWeight=0;
        TTree *tIn=0;    // input tree to read

        // objects to read from the tree
        panda::Event event;

        // configuration read from output tree
        std::vector<int> ibetas;
        std::vector<int> Ns; 
        std::vector<int> orders;


    };
}

#endif

