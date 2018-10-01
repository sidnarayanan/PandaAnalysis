#include "../interface/BTagTools.h"

using namespace pa;
using namespace std;

BTagCorrs::BTagCorrs(TString dirPath, const Analysis& analysis, GeneralTree& gt_) : 
  readers(bN),
  gt(gt_)
{
    for (auto& r : readers)
      r.reset(nullptr); 

    if (analysis.btagSFs) {
        if (analysis.useDeepCSV) {
          if (analysis.year==2016)
            calib.reset(new BTagCalibration("DeepCSV", (dirPath+"csv/DeepCSV_Moriond17_B_H.csv").Data()));
          else if (analysis.year==2017)
            calib.reset(new BTagCalibration("DeepCSV", (dirPath+"csv/DeepCSV_94XSF_V2_B_F.csv").Data()));
        } else {
          if (analysis.year==2016)
            calib.reset(new BTagCalibration("csvv2", (dirPath+"moriond17/CSVv2_Moriond17_B_H.csv").Data()));
          else if (analysis.year==2017)
            calib.reset(new BTagCalibration("csvv2", (dirPath+"csv/CSVv2_94XSF_V2_B_F.csv").Data()));
        }
        readers[bJetL].reset(new BTagCalibrationReader(BTagEntry::OP_LOOSE,"central",{"up","down"}));
        readers[bJetL]->load(*(calib),BTagEntry::FLAV_B,"comb");
        readers[bJetL]->load(*(calib),BTagEntry::FLAV_C,"comb");
        readers[bJetL]->load(*(calib),BTagEntry::FLAV_UDSG,"incl");
        
        if (analysis.year==2017) 
          sj_calib.reset(new BTagCalibration("DeepCSV",(dirPath+"csv/subjet_DeepCSV_94XSF_V2_B_F.csv").Data()));
        else 
          sj_calib.reset(new BTagCalibration("csvv2",(dirPath+"moriond17/subjet_CSVv2_Moriond17_B_H.csv").Data()));
        readers[bSubJetL].reset(new BTagCalibrationReader(BTagEntry::OP_LOOSE,"central",{"up","down"}));
        readers[bSubJetL]->load(*(sj_calib),BTagEntry::FLAV_B,"lt");
        readers[bSubJetL]->load(*(sj_calib),BTagEntry::FLAV_C,"lt");
        readers[bSubJetL]->load(*(sj_calib),BTagEntry::FLAV_UDSG,"incl");

        readers[bJetM].reset(new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"}));
        readers[bJetM]->load(*(calib),BTagEntry::FLAV_B,"comb");
        readers[bJetM]->load(*(calib),BTagEntry::FLAV_C,"comb");
        readers[bJetM]->load(*(calib),BTagEntry::FLAV_UDSG,"incl");
        
        readers[bSubJetM].reset(new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"}));
        readers[bSubJetM]->load(*(sj_calib),BTagEntry::FLAV_B,"lt");
        readers[bSubJetM]->load(*(sj_calib),BTagEntry::FLAV_C,"lt");
        readers[bSubJetM]->load(*(sj_calib),BTagEntry::FLAV_UDSG,"incl");
    }

    if (analysis.btagWeights) {
      // Build a vector of syst names for instantiating the BTagCalibrationReader
      std::vector<std::string> systNames;
      systNames.reserve(GeneralTree::nCsvShifts);
      for (unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
        GeneralTree::csvShift shift = gt.csvShifts[iShift];
        if (shift==GeneralTree::csvCent) 
          continue;
        systNames.push_back(GeneralTree::csvShiftName(shift).Data());
      }    
      reshaper.reset(new BTagCalibrationReader(BTagEntry::OP_RESHAPING,
                                               GeneralTree::csvShiftName(GeneralTree::csvCent).Data(), 
                                               systNames));

      if (analysis.year==2016) {
        if (analysis.useCMVA) {
          reshaper_calib.reset(new BTagCalibration("cMVAv2", (dirPath+"moriond17/cMVAv2_Moriond17_B_H.csv").Data()));
        } else {
          if (analysis.useDeepCSV) 
            // Waiting for the Official 2016 SFs!!!
            reshaper_calib.reset(new BTagCalibration("DeepCSV", (dirPath+"csv/DeepCSV_94XSF_V2_B_F.csv").Data()));
          else
            reshaper_calib.reset(new BTagCalibration("csvv2", (dirPath+"moriond17/CSVv2_Moriond17_B_H.csv").Data()));
        }
      } else if (analysis.year==2017) {
        if (analysis.useCMVA) {
          logger.warning("BTagCorrs::BTagCorrs","CMVA is not supported in 2017!");
          reshaper_calib.reset(new BTagCalibration("cMVAv2", (dirPath+"moriond17/cMVAv2_Moriond17_B_H.csv").Data()));
        } else {
          if (analysis.useDeepCSV)
            reshaper_calib.reset(new BTagCalibration("DeepCSV", (dirPath+"csv/DeepCSV_94XSF_V2_B_F.csv").Data()));
          else
            reshaper_calib.reset(new BTagCalibration("csvv2", (dirPath+"csv/CSVv2_94XSF_V2_B_F.csv").Data()));
        }
      }
      reshaper->load(*(reshaper_calib), BTagEntry::FLAV_B, "iterativeFit");
      reshaper->load(*(reshaper_calib), BTagEntry::FLAV_C, "iterativeFit");
      reshaper->load(*(reshaper_calib), BTagEntry::FLAV_UDSG, "iterativeFit");
    }
}

void BTagCorrs::calcSF(BTagType bt, int flavor,
                       double eta, double pt, double eff, double uncFactor,
                       double &sf, double &sfUp, double &sfDown) 
{
  if (flavor==5) {
    sf     = readers[bt]->eval_auto_bounds("central",BTagEntry::FLAV_B,eta,pt);
    sfUp   = readers[bt]->eval_auto_bounds("up",BTagEntry::FLAV_B,eta,pt);
    sfDown = readers[bt]->eval_auto_bounds("down",BTagEntry::FLAV_B,eta,pt);
  } else if (flavor==4) {
    sf     = readers[bt]->eval_auto_bounds("central",BTagEntry::FLAV_C,eta,pt);
    sfUp   = readers[bt]->eval_auto_bounds("up",BTagEntry::FLAV_C,eta,pt);
    sfDown = readers[bt]->eval_auto_bounds("down",BTagEntry::FLAV_C,eta,pt);
  } else {
    sf     = readers[bt]->eval_auto_bounds("central",BTagEntry::FLAV_UDSG,eta,pt);
    sfUp   = readers[bt]->eval_auto_bounds("up",BTagEntry::FLAV_UDSG,eta,pt);
    sfDown = readers[bt]->eval_auto_bounds("down",BTagEntry::FLAV_UDSG,eta,pt);
  }

  sfUp = uncFactor*(sfUp-sf)+sf;
  sfDown = uncFactor*(sfDown-sf)+sf;
  return;
}

void BTagCorrs::evalSF(vector<btagcand> &cands, vector<double> &sfs,
                       GeneralTree::BTagShift shift,GeneralTree::BTagJet jettype, bool do2) 
{
  float sf0 = 1, sf1 = 1, sfGT0 = 1, sf2=1;
  float prob_mc0=1, prob_data0=1;
  float prob_mc1=0, prob_data1=0;
  unsigned int nC = cands.size();

  for (unsigned int iC=0; iC!=nC; ++iC) {
    double sf_i = sfs[iC];
    double eff_i = cands[iC].eff;
    prob_mc0 *= (1-eff_i);
    prob_data0 *= (1-sf_i*eff_i);
    float tmp_mc1=1, tmp_data1=1;
    for (unsigned int jC=0; jC!=nC; ++jC) {
      if (iC==jC) continue;
      double sf_j = sfs[jC];
      double eff_j = cands[jC].eff;
      tmp_mc1 *= (1-eff_j);
      tmp_data1 *= (1-eff_j*sf_j);
    }
    prob_mc1 += eff_i * tmp_mc1;
    prob_data1 += eff_i * sf_i * tmp_data1;
  }
  
  if (nC>0) {
    sf0 = prob_data0/prob_mc0;
    sf1 = prob_data1/prob_mc1;
    sfGT0 = (1-prob_data0)/(1-prob_mc0);
  }

  GeneralTree::BTagParams p;
  p.shift = shift;
  p.jet = jettype;
  p.tag=GeneralTree::b0; gt.sf_btags[p] = sf0;
  p.tag=GeneralTree::b1; gt.sf_btags[p] = sf1;
  p.tag=GeneralTree::bGT0; gt.sf_btags[p] = sfGT0;

  if (do2) {
    float prob_mc2=0, prob_data2=0;
    unsigned int nC = cands.size();

    for (unsigned int iC=0; iC!=nC; ++iC) {
      double sf_i = sfs[iC], eff_i = cands[iC].eff;
      for (unsigned int jC=iC+1; jC!=nC; ++jC) {
        double sf_j = sfs[jC], eff_j = cands[jC].eff;
        float tmp_mc2=1, tmp_data2=1;
        for (unsigned int kC=0; kC!=nC; ++kC) {
          if (kC==iC || kC==jC) continue;
          double sf_k = sfs[kC], eff_k = cands[kC].eff;
          tmp_mc2 *= (1-eff_k);
          tmp_data2 *= (1-eff_k*sf_k);
        }
        prob_mc2 += eff_i * eff_j * tmp_mc2;
        prob_data2 += eff_i * sf_i * eff_j * sf_j * tmp_data2;
      }
    }

    if (nC>1) {
      sf2 = prob_data2/prob_mc2;
    }

    p.tag=GeneralTree::b2; gt.sf_btags[p] = sf2;
  }
}
