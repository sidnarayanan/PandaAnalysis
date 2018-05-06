#ifndef CONFIG
#define CONFIG

#include "TString.h"
#include "vector"
#include "map"

#include "AnalyzerUtilities.h"
#include "Common.h"

// fastjet
#include "fastjet/contrib/MeasureDefinition.hh"

// JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

namespace pa {

  class Config {
  public:
    Config(const Analysis& a_, int DEBUG_ = 0) : 
    DEBUG(DEBUG_),
    analysis(a_),
    tr("PandaAnalyzer::Run", DEBUG+1),
    isData(analysis.isData)
    { }


    const int DEBUG;
    const Analysis& analysis;
    TimeReporter tr;

    float minJetPt{30};
    float minBJetPt{30};
    int maxshiftJES{1};
    float minGenBosonPt{100};
    float maxGenBosonPt{1000};
    float minGenFatJetPt{450};
    float minSoftTrackPt{0.3};

    bool isData;        // to do gen matching, etc

    int NPFPROPS{9}, NSVPROPS{13}, NMAXPF{100}, NMAXSV{10}, NGENPROPS{8};
    float FATJETMATCHDR2{2.25};
    int NJETSAVED{2};


    // configuration read from output tree
    std::vector<int> ibetas;
    std::vector<int> Ns; 
    std::vector<int> orders;

    TString auxFilePath{""}; 
  };

  class Utils {
  public:
    Utils() { }
    ~Utils(); 

    double getCorr(CorrectionType ct, double x, double y=0);
    double getError(CorrectionType ct, double x, double y=0);
    void openCorr(CorrectionType, TString, TString, int);

    std::map<int, int> pdgToQ; 
    
    // fastjet reclustering
    fastjet::AreaDefinition    *areaDef         {nullptr};
    fastjet::GhostedAreaSpec   *activeArea      {nullptr};

    EraHandler *eras{nullptr{; //!< determining data-taking era, to be used for era-dependent JEC

    auto   fCorrs = std::vector<TFile*>  (cN,nullptr);
    auto  f1Corrs = std::vector<TF1Corr*>  (cN,nullptr);
    auto  h1Corrs = std::vector<TH1Corr*>  (cN,nullptr);
    auto  h2Corrs = std::vector<TH2Corr*>  (cN,nullptr);

    TFile   *fMSDcorr       {nullptr};
    TF1     *puppisd_corrGEN    {nullptr};
    TF1     *puppisd_corrRECO_cen {nullptr};
    TF1     *puppisd_corrRECO_for {nullptr};

    BTagCorrs *btag{nullptr};

  };

}

#endif
