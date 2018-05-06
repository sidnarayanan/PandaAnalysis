#ifndef CONFIG
#define CONFIG

#include "TString.h"
#include "vector"
#include "map"

#include "AnalyzerUtilities.h"
#include "Common.h"
#include "BTagTools.h"

// fastjet
#include "fastjet/contrib/MeasureDefinition.hh"
#include "fastjet/contrib/SoftDrop.hh"

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
    tr("PandaAnalyzer", DEBUG+1),
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
    Utils() : 
      fCorrs(cN,nullptr),
      f1Corrs(cN,nullptr),
      h1Corrs(cN,nullptr),
      h2Corrs(cN,nullptr)   { }
    ~Utils(); 

    double getCorr(CorrectionType ct, double x, double y=0);
    double getError(CorrectionType ct, double x, double y=0);
    void openCorr(CorrectionType, TString, TString, int);

    std::map<int, int> pdgToQ; 
    
    // fastjet reclustering
    fastjet::AreaDefinition    *areaDef         {nullptr};
    fastjet::GhostedAreaSpec   *activeArea      {nullptr};
    fastjet::contrib::SoftDrop   *softDrop      {nullptr};

    EraHandler *eras{nullptr}; //!< determining data-taking era, to be used for era-dependent JEC

    std::vector<TFile*> fCorrs;
    std::vector<TF1Corr*> f1Corrs;
    std::vector<THCorr1*> h1Corrs;
    std::vector<THCorr2*> h2Corrs;

    TFile   *fMSDcorr       {nullptr};
    TF1     *puppisd_corrGEN    {nullptr};
    TF1     *puppisd_corrRECO_cen {nullptr};
    TF1     *puppisd_corrRECO_for {nullptr};

    BTagCorrs *btag{nullptr};

  };

}

#endif
