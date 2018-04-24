#ifndef PandaAnalyzer_h
#define PandaAnalyzer_h

// STL
#include "vector"
#include <unordered_set>
#include "map"
#include <string>
#include <cmath>

// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TRandom3.h>
#include <TLorentzVector.h>

#include "AnalyzerUtilities.h"
#include "GeneralTree.h"
#include "Selection.h"

// fastjet
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/MeasureDefinition.hh"

// btag
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
//#include "BTagCalibrationStandalone.h"

// JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// Utils
#include "PandaAnalysis/Utilities/interface/RoccoR.h"
#include "PandaAnalysis/Utilities/interface/CSVHelper.h"
#include "PandaAnalysis/Utilities/interface/EnergyCorrelations.h"

// TMVA
#include "TMVA/Reader.h"

// macros
#define JESLOOP for (int shift = 0; shift != cfg.maxshiftJES; ++shift)

inline int jes2i(shiftjes i) { return static_cast<int>(i); }
inline shiftjes i2jes(int i) { return static_cast<shiftjes>(i); }

/////////////////////////////////////////////////////////////////////////////
// PandaAnalyzer definition
class PandaAnalyzer {
public :

    PandaAnalyzer(int debug_=0);
    ~PandaAnalyzer();
    int Init(TTree *tree, TH1D *hweights, TTree *weightNames=0);
    void SetOutputFile(TString fOutName);
    void ResetBranches();
    void Run();
    void Terminate();
    void SetDataDir(const char *s);
    void AddPresel(Selection *s) { selections.push_back(s); }
    void AddGoodLumiRange(int run, int l0, int l1);

    // public configuration
    void SetAnalysis(Analysis *a) { analysis = a; }

private:
    //////////////////////////////////////////////////////////////////////////////////////

    bool PassGoodLumis(int run, int lumi);
    bool PassPresel(Selection::Stage stage);
    void OpenCorrection(CorrectionType,TString,TString,int);
    double GetCorr(CorrectionType ct,double x, double y=0);
    double GetError(CorrectionType ct,double x, double y=0);
    void RegisterTriggers(); 
    
    // miscellaneous helper functions
    float getMSDCorr(float, float); 
    void jetPartonFlavor(const panda::Jet&, int& flavor, float& genpt);
    void jetClusteredFlavor(const panda::Jet&, int& flavor, float& genpt);
    JetWrapper shiftJet(const panda::Jet&, shiftjes);
    double weightEWKCorr(float pt, int type);
    double weightZHEWKCorr(float baseCorr);

    // these are functions used for analysis-specific tasks inside Run.
    // ideally the return type is void (e.g. they are stateful functions),
    void CalcBJetSFs(BTagType bt, int flavor, double eta, double pt, 
                     double eff, double uncFactor, double &sf, double &sfUp, double &sfDown);
    void ComplicatedLeptons();
    void ComplicatedPhotons();
    void EvalBTagSF(std::vector<btagcand> &cands, std::vector<double> &sfs,
                    GeneralTree::BTagShift shift,GeneralTree::BTagJet jettype, bool do2=false);
    void IncrementAuxFile(bool close=false);
    void IncrementGenAuxFile(bool close=false);
    void FatjetBasics();
    void FatjetMatching();
    void FatjetPartons();
    void FatjetRecluster();
    void FillPFTree();
    void GenFatJet();
    void GenJetsNu();
    void GenStudyEWK();
    void GetMETSignificance(); 
    void HeavyFlavorCounting();
    void IsoJet(JetWrapper&, JESHandler&);
    void JetBRegressionInfo(const panda::Jet&, int, int);
    void InclusiveLeptons();
    void JetBasics();
    void JetBtagSFs();
    void JetCMVAWeights();
    void JetHbbBasics(const panda::Jet&);
    void JetHbbReco(int);
    void JetHbbSoftActivity();
    void JetVBFBasics(const panda::Jet&);
    void JetVBFSystem(int);
    void JetVaryJES();
    void LeptonSFs();
    void LHEInfo();
    bool PFChargedPhotonMatch(const panda::Photon&);
    void PhotonSFs();
    void QCDUncs();
    void Recoil();
    void SaveGenLeptons();
    void SetupJES();
    void SignalInfo();
    void SignalReweights();
    void SimpleLeptons();
    void SimplePhotons();
    void Taus();
    void TopPTReweight();
    void TriggerEffs();
    void VJetsReweight();

    // templated functions
    // T = gen particle type
    template <typename T> void CountGenPartons(std::unordered_set<const T*>&);
    template <typename T> void FillGenTree();
    template <typename T> void RemoveGenDups(const panda::Collection<T>&);
    // T = gen jet type
    template <typename T> void MatchGenJets(T& genJets);

    //////////////////////////////////////////////////////////////////////////////////////

    int DEBUG = 0; //!< debug verbosity level
    Analysis *analysis = 0; //!< configure what to run

    //////////////////////////////////////////////////////////////////////////////////////

    // stuff for matching objects
    std::map<panda::GenParticle const*,float> genObjects; 
        //!< particles we want to match the jets to, and the 'size' of the daughters
    panda::GenParticle const* MatchToGen(double eta, double phi, double r2, int pdgid=0);   
        //!< private function to match a jet; returns NULL if not found
    std::map<int,std::vector<LumiRange>> goodLumis;
    std::vector<panda::Particle*> matchPhos, matchEles, matchLeps;

    //////////////////////////////////////////////////////////////////////////////////////

    std::vector<Selection*> selections;

    //////////////////////////////////////////////////////////////////////////////////////

    // files and histograms containing weights
    std::vector<TFile*>   fCorrs  = std::vector<TFile*>  (cN,0); //!< files containing corrections
    std::vector<THCorr1*> h1Corrs = std::vector<THCorr1*>(cN,0); //!< histograms for binned corrections
    std::vector<THCorr2*> h2Corrs = std::vector<THCorr2*>(cN,0); //!< histograms for binned corrections
    std::vector<TF1Corr*> f1Corrs = std::vector<TF1Corr*>(cN,0); //!< TF1s for continuous corrections
    TFile *MSDcorr=0;
    TF1* puppisd_corrGEN=0;
    TF1* puppisd_corrRECO_cen=0;
    TF1* puppisd_corrRECO_for=0;
    RoccoR *rochesterCorrection=0;
    CSVHelper *csvReweighter=0, *cmvaReweighter=0;

    //////////////////////////////////////////////////////////////////////////////////////

    // IO for the analyzer
    TString fOutPath;
    TFile *fOut=0;     // output file is owned by PandaAnalyzer
    TTree *tOut=0;
    GeneralTree *gt=0; // essentially a wrapper around tOut
    TString auxFilePath="";
    unsigned auxCounter=0;
    TFile *fAux=0; // auxillary file
    TTree *tAux=0;
    TH1D *hDTotalMCWeight=0;
    TTree *tIn=0;    // input tree to read
    unsigned int preselBits=0;
    panda::EventAnalysis event;

    //////////////////////////////////////////////////////////////////////////////////////

    // configuration read from output tree
    std::vector<int> ibetas;
    std::vector<int> Ns; 
    std::vector<int> orders;

    //////////////////////////////////////////////////////////////////////////////////////

    // stuff that gets passed between modules
    // NB: ensure that any global vectors/maps that are per-event
    // are reset properly in ResetBranches(), or you can really
    // mess up behavior
    std::vector<TString> wIDs;
    std::vector<TriggerHandler> triggerHandlers = std::vector<TriggerHandler>(kNTrig);

    std::vector<panda::Lepton*> looseLeps, tightLeps, inclusiveLeps;
    std::vector<panda::Photon*> loosePhos;
    int looseLep1PdgId, looseLep2PdgId, looseLep3PdgId, looseLep4PdgId;

    std::vector<const panda::Particle*> validGenP;

    const panda::FatJet *fj1 = 0;
    panda::FatJetCollection *fatjets = 0;

    float minJetPt=30;
    float minBJetPt=30;
    panda::JetCollection *ak4jets = 0;
    std::vector<JESHandler> jesShifts = std::vector<JESHandler>(jes2i(shiftjes::N)); 
    int maxshift = 1; // max JES shift

    std::vector<panda::GenJet> genJetsNu;
    float genBosonPtMin, genBosonPtMax;

    float *bjetreg_vars = 0;

    std::vector<std::vector<float>> pfInfo;
    std::vector<std::vector<float>> svInfo; 
    float fjmsd, fjpt, fjrawpt, fjeta, fjphi;
    int NPFPROPS = 9, NSVPROPS = 13;
    int NMAXPF = 100, NMAXSV = 10;

    GenJetInfo genJetInfo;
    int NGENPROPS = 8; 
    float minGenFatJetPt = 450;
    
    float minSoftTrackPt = 0.3; // 300 MeV

};


/** templated functions **/

#include "TemplatedPandaAnalyzer.h"

#endif

