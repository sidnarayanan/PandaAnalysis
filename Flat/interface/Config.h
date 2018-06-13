#ifndef CONFIG
#define CONFIG

#include "TString.h"
#include "vector"
#include "map"
#include "memory.h"

#include "AnalyzerUtilities.h"
#include "PandaAnalysis/Flat/interface/Common.h"
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
    { 
      tr.Start();
    }


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
    Utils() : fCorrs(cN), f1Corrs(cN), h1Corrs(cN), h2Corrs(cN) { 
      for (unsigned i = 0; i != cN; ++i) {
        fCorrs[i].reset(nullptr);
        f1Corrs[i].reset(nullptr);
        h1Corrs[i].reset(nullptr);
        h2Corrs[i].reset(nullptr);
      }
    }
    ~Utils(); 

    double getCorr(CorrectionType ct, double x, double y=0);
    double getError(CorrectionType ct, double x, double y=0);
    void openCorr(CorrectionType, TString, TString, int);

    std::map<int, int> pdgToQ; 
    
    // fastjet reclustering
    std::unique_ptr<fastjet::AreaDefinition> areaDef{nullptr};
    std::unique_ptr<fastjet::GhostedAreaSpec> activeArea{nullptr};
    std::unique_ptr<fastjet::contrib::SoftDrop> softDrop{nullptr};

    std::unique_ptr<EraHandler> eras{nullptr}; //!< determining data-taking era, to be used for era-dependent JEC

    std::vector<std::unique_ptr<TFile>> fCorrs;
    std::vector<std::unique_ptr<TF1Corr>> f1Corrs;
    std::vector<std::unique_ptr<THCorr1>> h1Corrs;
    std::vector<std::unique_ptr<THCorr2>> h2Corrs;

    std::unique_ptr<TFile>fMSDcorr{nullptr};
    TF1     *puppisd_corrGEN    {nullptr};
    TF1     *puppisd_corrRECO_cen {nullptr};
    TF1     *puppisd_corrRECO_for {nullptr};

    std::unique_ptr<BTagCorrs>btag{nullptr};

  };

}

#endif
