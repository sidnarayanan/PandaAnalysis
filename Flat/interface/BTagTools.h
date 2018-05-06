#ifndef BTAGTOOLS
#define BTAGTOOLS

// btag
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "PandaAnalysis/Utilities/interface/CSVHelper.h"

namespace pa {
  
  enum BTagType {
    bJetL=0,
    bSubJetL,
    bJetM,
    bSubJetM,
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

  class BTagCorrs {
    public:
      BTagCorrs(TString dirPath, const Analysis& a, GeneralTree& gt_);
      ~BTagCorrs() { 
        for (auto* r : readers)
          delete r;
        delete calib;
        delete sj_calib;
        delete csvReweighter;
        delete cmvaReweighter;
      }
    
      void evalSF(std::vector<btagcand> &cands, std::vector<double> &sfs,
                  GeneralTree::BTagShift shift, GeneralTree::BTagJet jettype, bool do2=false);
      void calcSF(BTagType bt, int flavor, double eta, double pt, 
                  double eff, double uncFactor, double &sf, double &sfUp, double &sfDown);
      std::vector<BTagCalibrationReader*> readers = std::vector<BTagCalibrationReader*>(bN,nullptr);
      CSVHelper *csvReweighter{nullptr}, *cmvaReweighter{nullptr};
    private:
      BTagCalibration *calib{nullptr}, *sj_calib{nullptr};
      GeneralTree& gt; 
  };

}

#endif
