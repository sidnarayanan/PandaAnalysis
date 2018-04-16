#ifndef CSVHelper_h__
#define CSVHelper_h__

#include "TFile.h"
#include "TH1D.h"

// Adapted from ttH guys on Nov 21, 2017 ~DGH
// https://github.com/cms-ttH/MiniAOD/tree/master/MiniAODHelper

class CSVHelper
{
  public:
    // nHFptBins specifies how many of these pt bins are used:
    // (jetPt >= 19.99 && jetPt < 30), (jetPt >= 30 && jetPt < 40), (jetPt >= 40 && jetPt < 60), 
    // (jetPt >= 60 && jetPt < 100), (jetPt >= 100 && jetPt < 160), (jetPt >= 160 && jetPt < 10000).
    // If nHFptBins < 6, the last on is inclusive (eg jetPt >=100 && jetPt < 10000 for nHFptBins=5).
    // The SFs from data have 5 bins, the pseudo data scale factors 6 bins.
    CSVHelper(std::string hf="", std::string lf="", int nHFptBins=6);

    double getCSVWeight(std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs,
                       std::vector<int> jetFlavors, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF);

  private:
    void fillCSVHistos(TFile *fileHF, TFile *fileLF);

    // CSV reweighting
    TH1D *h_csv_wgt_hf[9][6];
    TH1D *c_csv_wgt_hf[9][6];
    TH1D *h_csv_wgt_lf[9][4][3];
    const int nHFptBins;
};

#endif
