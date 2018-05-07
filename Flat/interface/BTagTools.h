#ifndef BTAGTOOLS
#define BTAGTOOLS

// btag
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "PandaAnalysis/Utilities/interface/CSVHelper.h"

#include "Common.h"
#include "GeneralTree.h"

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
        delete reshaper_calib;
      }
    
      void evalSF(std::vector<btagcand> &cands, std::vector<double> &sfs,
                  GeneralTree::BTagShift shift, GeneralTree::BTagJet jettype, bool do2=false);
      void calcSF(BTagType bt, int flavor, double eta, double pt, 
                  double eff, double uncFactor, double &sf, double &sfUp, double &sfDown);
      std::vector<BTagCalibrationReader*> readers = std::vector<BTagCalibrationReader*>(bN,nullptr);
      BTagCalibrationReader *reshaper{nullptr};
      
      std::map<GeneralTree::csvShift, std::string> officialShiftNames = {
        {GeneralTree::csvCent         , "central"      },
        {GeneralTree::csvJESup        , "up_jes"       },
        {GeneralTree::csvJESdown      , "down_jes"     },
        {GeneralTree::csvLFup         , "up_lf"        },
        {GeneralTree::csvLFdown       , "down_lf"      },
        {GeneralTree::csvHFup         , "up_hf",       },
        {GeneralTree::csvHFdown       , "down_hf",     },
        {GeneralTree::csvHFStats1up   , "up_hfstats1"  },
        {GeneralTree::csvHFStats1down , "down_hfstats1"},
        {GeneralTree::csvHFStats2up   , "up_hfstats2"  },
        {GeneralTree::csvHFStats2down , "down_hfstats2"},
        {GeneralTree::csvLFStats1up   , "up_lfstats1"  },
        {GeneralTree::csvLFStats1down , "down_lfstats1"},
        {GeneralTree::csvLFStats2up   , "up_lfstats2"  },
        {GeneralTree::csvLFStats2down , "down_lfstats2"},
        {GeneralTree::csvCErr1up      , "up_cferr1"    },
        {GeneralTree::csvCErr1down    , "down_cferr1"  },
        {GeneralTree::csvCErr2up      , "up_cferr2"    },
        {GeneralTree::csvCErr2down    , "down_cferr2"  }
      };
    private:
      BTagCalibration *calib{nullptr}, *sj_calib{nullptr}, *reshaper_calib{nullptr};
      GeneralTree& gt;

  };

}

#endif
