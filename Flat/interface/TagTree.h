#ifndef TagTREE_H
#define TagTREE_H 

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TString.h"
#include "genericTree.h"
#include <map>

class TagTree : public genericTree {
    public:
      // public objects
      struct ECFParams {
        int ibeta;
        int N;
        int order;
        bool operator==(const ECFParams &o) const {
          return ibeta==o.ibeta && N==o.N && order==o.order;
        }
        bool operator<(const ECFParams &o) const {
          return ( N<o.N ||
                   (N==o.N && order<o.order) ||
                   (N==o.N && order==o.order && ibeta<o.ibeta) );
        }
        bool operator>(const ECFParams &o) const {
          return ! operator<(o);
        }
      };

        
    private:
        std::vector<double> betas = {0.5, 1.0, 2.0, 4.0};
        std::vector<int> ibetas = {0,1,2,3};
        std::vector<int> Ns = {1,2,3,4}; 
        std::vector<int> orders = {1,2,3};
        std::vector<ECFParams> ecfParams;

        TString makeECFString(ECFParams p) {
            return TString::Format("ECFN_%i_%i_%.2i",p.order,p.N,int(10*betas.at(p.ibeta)));
        }

    public:
      TagTree();
      ~TagTree();
      void WriteTree(TTree *t);
      void Fill() { treePtr->Fill(); }
      void Reset();

      std::vector<double> get_betas() const { return betas; }
      std::vector<int> get_ibetas() const { return ibetas; }
      std::vector<int> get_Ns() const { return Ns; }
      std::vector<int> get_orders() const { return orders; }
        

//STARTCUSTOMDEF
      std::map<ECFParams,float> fj1ECFNs;
//ENDCUSTOMDEF
    int runNumber = -1;
    int lumiNumber = -1;
    ULong64_t eventNumber = -1;
    int npv = -1;
    int pu = -1;
    float mcWeight = -1;
    float filter_maxRecoil = -1;
    float filter_whichRecoil = -1;
    float sf_ewkV = -1;
    float sf_qcdV = -1;
    float sf_ewkV2j = -1;
    float sf_qcdV2j = -1;
    float sf_qcdTT = -1;
    float sf_pu = -1;
    float sf_npv = -1;
    float sf_tt = -1;
    float sf_phoPurity = -1;
    float pfmet = -1;
    float puppimet = -1;
    float partonPt = -1;
    float partonEta = -1;
    int partonIsHad = -1;
    float partonSize = -1;
    int partonIsReco = -1;
    int partonPdgId = -1;
    int nFatJet = -1;
    float fj1Tau32 = -1;
    float fj1Tau21 = -1;
    float fj1Tau32SD = -1;
    float fj1Tau21SD = -1;
    float fj1MSD = -1;
    float fj1Pt = -1;
    float fj1Phi = -1;
    float fj1Eta = -1;
    float fj1M = -1;
    float fj1MaxCSV = -1;
    float fj1SubMaxCSV = -1;
    float fj1MinCSV = -1;
    float fj1DoubleCSV = -1;
    int fj1IsTight = -1;
    int fj1IsLoose = -1;
    float fj1RawPt = -1;
    int fj1NHF = -1;
    float fj1HTTMass = -1;
    float fj1HTTFRec = -1;
    int fj1IsClean = -1;
};

#endif
