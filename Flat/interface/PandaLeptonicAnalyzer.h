#ifndef PandaLeptonicAnalyzer_h
#define PandaLeptonicAnalyzer_h

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
#include "GeneralLeptonicTree.h"

// btag
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
//#include "BTagCalibrationStandalone.h"

// JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

/////////////////////////////////////////////////////////////////////////////
// some misc definitions


/////////////////////////////////////////////////////////////////////////////
// PandaLeptonicAnalyzer definition
class PandaLeptonicAnalyzer {
public :
    // configuration enums
    enum PreselectionBit {
     kTriggers   =(1<<0),
     kLepton     =(1<<1),
     kLeptonFake =(1<<2)
    };

    enum SelectionBit {
     kLoose       =(1<<0),
     kFake        =(1<<1),
     kMedium      =(1<<2),
     kTight       =(1<<3),
     kDxyz        =(1<<4),
     kTrigger     =(1<<5),
     kFakeTrigger =(1<<6)
    };

    enum ProcessType { 
        kNone,
        kZPtCut,
        kW,
        kA,
        kZEWK,
        kWEWK,
        kTT,
        kTop, // used for non-ttbar top
        kV, // used for non V+jets W or Z
        kH,
        kSignal,
    };

    enum TriggerBits {
        kMETTrig       =(1<<0),
        kSinglePhoTrig =(1<<1),
        kMuEGTrig      =(1<<2),
        kMuMuTrig      =(1<<3),
        kMuTrig        =(1<<4),
        kEGEGTrig      =(1<<5),
        kEGTrig        =(1<<6),
        kMuTagTrig     =(1<<7),
        kEGTagTrig     =(1<<8),
        kMuFakeTrig    =(1<<9),
        kEGFakeTrig    =(1<<10)
    };

    //////////////////////////////////////////////////////////////////////////////////////

    PandaLeptonicAnalyzer(int debug_=0);
    ~PandaLeptonicAnalyzer();
    int Init(TTree *tree, TH1D *hweights, TTree *weightNames=0);
    void SetOutputFile(TString fOutName);
    void ResetBranches();
    void Run();
    void Terminate();
    void SetDataDir(const char *s2);
    double phi_star_eta(TLorentzVector lep1, TLorentzVector lep2, int pdgId1);
    void SetPreselectionBit(PreselectionBit b,bool on=true) {
        if (on) 
            preselBits |= b;
        else 
            preselBits &= ~b;
    }
    void AddGoodLumiRange(int run, int l0, int l1);

    // public configuration
    void SetFlag(TString flag, bool b=true) { flags[flag]=b; }
    bool isData=false;                                                 // to do gen matching, etc
    int firstEvent=-1;
    int lastEvent=-1;                                                    // max events to process; -1=>all
    ProcessType processType=kNone;                         // determine what to do the jet matching to

private:
    enum CorrectionType { //!< enum listing relevant corrections applied to MC
        cPU=0,         //!< true pu weight
        cPUUp,         //!< true pu weight Up
        cPUDown,       //!< true pu weight Down
	ZHEwkCorr,     //!< ZH Ewk Corr weight  
	ZHEwkCorrUp,   //!< ZH Ewk Corr weight Up  
	ZHEwkCorrDown, //!< ZH Ewk Corr weight Down  
	cLooseMuonId,
	cMediumMuonId,
	cTightMuonId,
	cLooseMuonIso, 
	cMediumMuonIso,
	cTightMuonIso,
	cTrackingMuon,
	cLooseElectronId, 
	cMediumElectronId,
	cTightElectronId,
	cTrackingElectron,
        cN
    };

    enum BTagType {
        bJetL=0,
        bJetM,
        bN
    };

    class btagcand {
        public:
            btagcand(unsigned int i, int f,double e,double cent,double up,double down) {
                idx = i;
                flav = f;
                eff = e;
                sf = cent;
                sfup = up;
                sfdown = down;
            }
            ~btagcand() { }
            int flav, idx;
            double eff, sf, sfup, sfdown;
    };

    const int nBinEta = 28;
    const int nBinPt = 36;
    const int nBinHighPt = 5;
    const int nBinRap = 12;
    const int nBinPhiStar = 34;
    const int nBinPtRap0 = 34;
    const int nBinPtRap1 = 34;
    const int nBinPtRap2 = 34;
    const int nBinPtRap3 = 34;
    const int nBinPtRap4 = 34;
    const int nBinWWSS = 5;
    const int nBinWW = 9;
    const int nBinWWN0JET = 3;
    bool PassGoodLumis(int run, int lumi);
    bool PassPreselection();
    void CalcBJetSFs(BTagType bt, int flavor, double eta, double pt, 
                         double eff, double uncFactor, double &sf, double &sfUp, double &sfDown);
    void EvalBTagSF(std::vector<btagcand> &cands, std::vector<double> &sfs,
                    GeneralLeptonicTree::BTagShift shift,GeneralLeptonicTree::BTagJet jettype, bool do2=false);
    void OpenCorrection(CorrectionType,TString,TString,int);
    double GetCorr(CorrectionType ct,double x, double y=0);
    double GetError(CorrectionType ct,double x, double y=0);
    void RegisterTrigger(TString path, std::vector<unsigned> &idxs); 

    int DEBUG = 0; //!< debug verbosity level
    std::map<TString,bool> flags;

    std::map<panda::GenParticle const*,float> genObjects;                 //!< particles we want to match the jets to, and the 'size' of the daughters
    panda::GenParticle const* MatchToGen(double eta, double phi, double r2, int pdgid=0);        //!< private function to match a jet; returns NULL if not found
    std::map<int,std::vector<LumiRange>> goodLumis;
    std::vector<panda::Particle*> matchPhos, matchEles, matchLeps;
    
    // CMSSW-provided utilities

    BTagCalibration *btagCalib=0;

    std::vector<BTagCalibrationReader*> btagReaders = std::vector<BTagCalibrationReader*>(bN,0); //!< maps BTagType to a reader 
    
    std::map<TString,JetCorrectionUncertainty*> ak8UncReader; //!< calculate JES unc on the fly
    JERReader *ak8JERReader; //!< fatjet jet energy resolution reader
    std::map<TString,JetCorrectionUncertainty*> ak4UncReader; //!< calculate JES unc on the fly
    std::map<TString,FactorizedJetCorrector*> ak4ScaleReader; //!< calculate JES on the fly
    JERReader *ak4JERReader; //!< fatjet jet energy resolution reader
    EraHandler eras = EraHandler(2016); //!< determining data-taking era, to be used for era-dependent JEC

    // files and histograms containing weights
    std::vector<TFile*> fCorrs = std::vector<TFile*>(cN,0); //!< files containing corrections
    std::vector<THCorr1*> h1Corrs = std::vector<THCorr1*>(cN,0); //!< histograms for binned corrections
    std::vector<THCorr2*> h2Corrs = std::vector<THCorr2*>(cN,0); //!< histograms for binned corrections

    // IO for the analyzer
    TFile *fOut;     // output file is owned by PandaLeptonicAnalyzer
    TTree *tOut;
    GeneralLeptonicTree *gt; // essentially a wrapper around tOut
    TH1F *hDTotalMCWeight=0;
/*
    TH1D *hDReco_Eta;
    TH1D *hDRecoMuon_P[28];
    TH1D *hDRecoMuon_F[28];
    TH1D *hDRecoTrack_P[28];
    TH1D *hDRecoTrack_F[28];
    TH1D *hDRecoMuonIso_P[28];
    TH1D *hDRecoMuonIso_F[28];
    TH1D *hDRecoTrackIso_P[28];
    TH1D *hDRecoTrackIso_F[28];
    TH1D *hDGenToMuon[8];
*/
    TH1D *hDDilPtMM;     TH1D *hDDilPtMM_PDF;	  TH1D *hDDilPtMM_QCD;	   TH1D *hDDilPtMM_QCDPart[6];
    TH1D *hDDilPtEE;	 TH1D *hDDilPtEE_PDF;	  TH1D *hDDilPtEE_QCD;	   TH1D *hDDilPtEE_QCDPart[6];
    TH1D *hDDilHighPtMM; TH1D *hDDilHighPtMM_PDF; TH1D *hDDilHighPtMM_QCD; TH1D *hDDilHighPtMM_QCDPart[6];
    TH1D *hDDilHighPtEE; TH1D *hDDilHighPtEE_PDF; TH1D *hDDilHighPtEE_QCD; TH1D *hDDilHighPtEE_QCDPart[6];
    TH1D *hDDilRapMM;	 TH1D *hDDilRapMM_PDF;	  TH1D *hDDilRapMM_QCD;    TH1D *hDDilRapMM_QCDPart[6];
    TH1D *hDDilRapEE;	 TH1D *hDDilRapEE_PDF;	  TH1D *hDDilRapEE_QCD;	   TH1D *hDDilRapEE_QCDPart[6];
    TH1D *hDDilPhiStarMM;TH1D *hDDilPhiStarMM_PDF;TH1D *hDDilPhiStarMM_QCD;TH1D *hDDilPhiStarMM_QCDPart[6];
    TH1D *hDDilPhiStarEE;TH1D *hDDilPhiStarEE_PDF;TH1D *hDDilPhiStarEE_QCD;TH1D *hDDilPhiStarEE_QCDPart[6];
    TH1D *hDDilPtRap0MM; TH1D *hDDilPtRap0MM_PDF; TH1D *hDDilPtRap0MM_QCD; TH1D *hDDilPtRap0MM_QCDPart[6];
    TH1D *hDDilPtRap0EE; TH1D *hDDilPtRap0EE_PDF; TH1D *hDDilPtRap0EE_QCD; TH1D *hDDilPtRap0EE_QCDPart[6];
    TH1D *hDDilPtRap1MM; TH1D *hDDilPtRap1MM_PDF; TH1D *hDDilPtRap1MM_QCD; TH1D *hDDilPtRap1MM_QCDPart[6];
    TH1D *hDDilPtRap1EE; TH1D *hDDilPtRap1EE_PDF; TH1D *hDDilPtRap1EE_QCD; TH1D *hDDilPtRap1EE_QCDPart[6];
    TH1D *hDDilPtRap2MM; TH1D *hDDilPtRap2MM_PDF; TH1D *hDDilPtRap2MM_QCD; TH1D *hDDilPtRap2MM_QCDPart[6];
    TH1D *hDDilPtRap2EE; TH1D *hDDilPtRap2EE_PDF; TH1D *hDDilPtRap2EE_QCD; TH1D *hDDilPtRap2EE_QCDPart[6];
    TH1D *hDDilPtRap3MM; TH1D *hDDilPtRap3MM_PDF; TH1D *hDDilPtRap3MM_QCD; TH1D *hDDilPtRap3MM_QCDPart[6];
    TH1D *hDDilPtRap3EE; TH1D *hDDilPtRap3EE_PDF; TH1D *hDDilPtRap3EE_QCD; TH1D *hDDilPtRap3EE_QCDPart[6];
    TH1D *hDDilPtRap4MM; TH1D *hDDilPtRap4MM_PDF; TH1D *hDDilPtRap4MM_QCD; TH1D *hDDilPtRap4MM_QCDPart[6];
    TH1D *hDDilPtRap4EE; TH1D *hDDilPtRap4EE_PDF; TH1D *hDDilPtRap4EE_QCD; TH1D *hDDilPtRap4EE_QCDPart[6];
    TH1D *hDWWSSMLL;
    TH1D *hDWWMLL;       TH1D *hDWWMLL_PDF;	  TH1D *hDWWMLL_QCD;	   TH1D *hDWWMLL_QCDPart[6];
    TH1D *hDWWPTL1;      TH1D *hDWWPTL1_PDF;	  TH1D *hDWWPTL1_QCD;	   TH1D *hDWWPTL1_QCDPart[6];
    TH1D *hDWWPTL2;      TH1D *hDWWPTL2_PDF;	  TH1D *hDWWPTL2_QCD;	   TH1D *hDWWPTL2_QCDPart[6];
    TH1D *hDWWDPHILL;    TH1D *hDWWDPHILL_PDF;	  TH1D *hDWWDPHILL_QCD;	   TH1D *hDWWDPHILL_QCDPart[6];
    TH1D *hDWWN0JET;     TH1D *hDWWN0JET_PDF;	  TH1D *hDWWN0JET_QCD;	   TH1D *hDWWN0JET_QCDPart[6];
    TH1D *hDWWNJET;      TH1D *hDWWNJET_PDF;	  TH1D *hDWWNJET_QCD;	   TH1D *hDWWNJET_QCDPart[6];
    TTree *tIn=0;    // input tree to read
    unsigned int preselBits=0;

    // objects to read from the tree
    panda::Event event;

    // any extra signal weights we want
    std::vector<TString> wIDs;

};

#endif
