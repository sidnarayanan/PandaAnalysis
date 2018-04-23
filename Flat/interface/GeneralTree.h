#ifndef GeneralTree_H
#define GeneralTree_H
// STARTCUSTOM HEADER

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TString.h"
#include "genericTree.h"
#include <map>

// ENDCUSTOM
#define NJET 20
#define NSUBJET 2
enum class shiftjes { 
  kNominal,
  kJESUp,
  kJESDown,
  N
}; 
#define NLEP 4
class GeneralTree : public genericTree {
  public:
    GeneralTree();
    ~GeneralTree();
    void WriteTree(TTree* t);
    void Fill() { treePtr->Fill(); }
    void Reset();
// STARTCUSTOM PUBLIC
    const std::vector<double>& get_betas() const { return betas; }
    const std::vector<int>& get_ibetas() const { return ibetas; }
    const std::vector<int>& get_Ns() const { return Ns; }
    const std::vector<int>& get_orders() const { return orders; }


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
    enum BTagShift {
      bCent=0,
      bBUp,
      bBDown,
      bMUp,
      bMDown,
      bNShift
    };
    enum BTagJet {
      bJet=0,
      bSubJet,
      bNJet
    };
    enum BTagTags {
      b0=0,
      b1,
      b2,
      bGT0,
      bNTags
    };
    struct BTagParams {
      BTagJet jet;
      BTagTags tag;
      BTagShift shift=(BTagShift)0;
      bool operator==(const BTagParams &o) const {
        return jet==o.jet && tag==o.tag && shift==o.shift;
      }
      bool operator<(const BTagParams &o) const {
        return ( jet<o.jet ||
                 (jet==o.jet && tag<o.tag) ||
                 (jet==o.jet && tag==o.tag && shift<o.shift) );
      }
      bool operator>(const BTagParams &o) const {
        return ! operator<(o);
      }
    };
    enum csvShift {
      csvCent=0,
      csvJESup=7,
      csvJESdown=8,
      csvLFup=9,
      csvLFdown=10,
      csvHFup=11,
      csvHFdown=12,
      csvHFStats1up=13,
      csvHFStats1down=14,
      csvHFStats2up=15,
      csvHFStats2down=16,
      csvLFStats1up=17,
      csvLFStats1down=18,
      csvLFStats2up=19,
      csvLFStats2down=20,
      csvCErr1up=21,
      csvCErr1down=22,
      csvCErr2up=23,
      csvCErr2down=24
    };
    // Array of the CSV/CMVA weight enums that can be looped over
    static const unsigned char nCsvShifts=19;
    csvShift csvShifts[nCsvShifts] = {
      csvCent,
      csvJESup,
      csvJESdown,
      csvLFup,
      csvLFdown,
      csvHFup,
      csvHFdown,
      csvHFStats1up,
      csvHFStats1down,
      csvHFStats2up,
      csvHFStats2down,
      csvLFStats1up,
      csvLFStats1down,
      csvLFStats2up,
      csvLFStats2down,
      csvCErr1up,
      csvCErr1down,
      csvCErr2up,
      csvCErr2down
    };
    virtual void SetAuxTree(TTree *t);

    // public config
    bool is_monohiggs=false, is_vbf=false, is_fatjet=true, 
         is_leptonic=false, is_photonic=false, is_monotop=true,
         is_hbb=false;
    bool btagWeights=false, useCMVA=false;
    std::map<ECFParams,float> fjECFNs;
    std::map<BTagParams,float> sf_btags;
    std::map<BTagParams,float> sf_alt_btags;
    std::map<TString,float> signal_weights;
    std::map<csvShift,float> sf_csvWeights;
// ENDCUSTOM
  private:
// STARTCUSTOM PRIVATE
    const std::vector<double> betas = {0.5, 1.0, 2.0, 4.0};
    const std::vector<int> ibetas = {0,1,2,3};
    const std::vector<int> Ns = {1,2,3,4};
    const std::vector<int> orders = {1,2,3};
    std::vector<ECFParams> ecfParams;
    std::vector<BTagParams> btagParams;
    TString makeECFString(ECFParams p) {
      return TString::Format("ECFN_%i_%i_%.2i",p.order,p.N,int(10*betas.at(p.ibeta)));
    }
    TString makeBTagSFString(BTagParams p) {
      TString s = "sf_";
      if (p.jet==bSubJet)
        s += "sj";
      s += "btag";
      switch (p.tag) {
        case b0:
          s += "0"; break;
        case b1:
          s += "1"; break;
        case b2:
          s += "2"; break;
        case bGT0:
          s += "GT0"; break;
        default:
          break;
      }
      if (p.shift==bCent)
        return s;
      switch (p.shift) {
        case bBUp:
          s += "BUp"; break;
        case bBDown:
          s += "BDown"; break;
        case bMUp:
          s += "MUp"; break;
        case bMDown:
          s += "MDown"; break;
        default: break;
      }
      return s;
    }
    TString makeCsvWeightString(csvShift shift, bool isCMVA=false) {
      TString s = isCMVA? "sf_cmvaWeight_" : "sf_csvWeight_";
      switch (shift) {
        case csvCent         :  s += "Cent"         ; break;
        case csvJESup        :  s += "JESup"        ; break;
        case csvJESdown      :  s += "JESdown"      ; break;
        case csvLFup         :  s += "LFup"         ; break;
        case csvLFdown       :  s += "LFdown"       ; break;
        case csvHFup         :  s += "HFup"         ; break;
        case csvHFdown       :  s += "HFdown"       ; break;
        case csvHFStats1up   :  s += "HFStats1up"   ; break;
        case csvHFStats1down :  s += "HFStats1down" ; break;
        case csvHFStats2up   :  s += "HFStats2up"   ; break;
        case csvHFStats2down :  s += "HFStats2down" ; break;
        case csvLFStats1up   :  s += "LFStats1up"   ; break;
        case csvLFStats1down :  s += "LFStats1down" ; break;
        case csvLFStats2up   :  s += "LFStats2up"   ; break;
        case csvLFStats2down :  s += "LFStats2down" ; break;
        case csvCErr1up      :  s += "CErr1up"      ; break;
        case csvCErr1down    :  s += "CErr1down"    ; break;
        case csvCErr2up      :  s += "CErr2up"      ; break;
        case csvCErr2down    :  s += "CErr2down"    ; break;
        default              :  s += "Unknown"      ; break;
      }
      return s;
    }
// ENDCUSTOM
  public:
  int runNumber;
  int lumiNumber;
  ULong64_t eventNumber;
  int isData;
  int npv;
  int pu;
  float mcWeight;
  int trigger;
  int metFilter;
  int egmFilter;
  float filter_maxRecoil;
  int filter_whichRecoil;
  int badECALFilter;
  float sf_ewkV;
  float sf_qcdV;
  float sf_ewkV2j;
  float sf_qcdV2j;
  float sf_qcdV_VBF;
  float sf_qcdV_VBF2l;
  float sf_qcdV_VBFTight;
  float sf_qcdV_VBF2lTight;
  float sf_qcdTT;
  float sf_lepID;
  float sf_lepIso;
  float sf_lepTrack;
  float sf_pho;
  float sf_eleTrig;
  float sf_muTrig;
  float sf_phoTrig;
  float sf_metTrig;
  float sf_metTrigZmm;
  float sf_metTrigVBF;
  float sf_metTrigZmmVBF;
  float sf_pu;
  float sf_npv;
  float sf_tt;
  float sf_phoPurity;
  float sumETRaw;
  float pfmetRaw;
  float pfmet[3];
  float pfmetphi[3];
  float pfmetnomu[3];
  float puppimet[3];
  float puppimetphi[3];
  float calomet;
  float calometphi;
  float pfcalobalance;
  float sumET;
  float trkmet;
  float trkmetphi;
  float pfmetsig;
  float puppimetsig;
  int whichRecoil;
  float puppiUWmag[3];
  float puppiUZmag[3];
  float puppiUAmag[3];
  float puppiUperp[3];
  float puppiUpara[3];
  float puppiUmag[3];
  float pfUWmag[3];
  float pfUZmag[3];
  float pfUAmag[3];
  float pfUperp[3];
  float pfUpara[3];
  float pfUmag[3];
  float puppiUWphi[3];
  float puppiUZphi[3];
  float puppiUAphi[3];
  float puppiUphi[3];
  float pfUWphi[3];
  float pfUZphi[3];
  float pfUAphi[3];
  float pfUphi[3];
  float dphipfmet[3];
  float dphipuppimet[3];
  float dphipuppiUW[3];
  float dphipuppiUZ[3];
  float dphipuppiUA[3];
  float dphipfUW[3];
  float dphipfUZ[3];
  float dphipfUA[3];
  float dphipuppiU[3];
  float dphipfU[3];
  float trueGenBosonPt;
  float genBosonPt;
  float genBosonEta;
  float genBosonMass;
  float genBosonPhi;
  float genMuonPt;
  float genMuonEta;
  float genElectronPt;
  float genElectronEta;
  float genTauPt;
  float genTauEta;
  float genJet1Pt;
  float genJet2Pt;
  float genJet1Eta;
  float genJet2Eta;
  float genMjj;
  float genTopPt;
  float genAntiTopPt;
  int nJet[3];
  int nJot[3];
  int nJotMax;
  int nIsoJet[3];
  float jetPt[2];
  float jetEta[2];
  float jetPhi[2];
  float jetGenPt[2];
  float jetCSV[2];
  int jetFlav[2];
  int jetIsTight[2];
  int jetIsIso[2];
  float jotPt[3][20];
  float jotE[3][20];
  float jotEta[3][20];
  float jotPhi[3][20];
  float jotCSV[3][20];
  int jotVBFID[3][20];
  float jotM[3][20];
  float jotCMVA[3][20];
  int jotIso[3][20];
  float jotQGL[3][20];
  float jotLep1Pt[3][20];
  float jotLep1PtRel[3][20];
  float jotLep1DeltaR[3][20];
  float jotTrk1Pt[3][20];
  float jotVtxPt[3][20];
  float jotVtxMass[3][20];
  float jotVtx3DVal[3][20];
  float jotVtx3DErr[3][20];
  int jotVtxNtrk[3][20];
  float jotEMF[3][20];
  float jotHF[3][20];
  int jotNLep[3][20];
  float jotGenPt[3][20];
  int jotFlav[3][20];
  int hbbjtidx[3][2];
  float hbbjtRegFac[3][2];
  float barrelJet1Pt;
  float barrelJet1Eta;
  float barrelHT;
  float barrelHTMiss;
  float barrelJet12Pt;
  float jot12Mass[3];
  float jot12DEta[3];
  float jot12DPhi[3];
  int jetNBtags[3];
  int jetNMBtags[3];
  int isojetNBtags[3];
  int nFatjet;
  float fjTau32;
  float fjTau21;
  float fjTau32SD;
  float fjTau21SD;
  float fjMSD;
  float fjRho;
  float fjRawRho;
  float fjRho2;
  float fjRawRho2;
  float fjMSD_corr;
  float fjPt;
  float fjPhi;
  float fjEta;
  float fjM;
  float fjMaxCSV;
  float fjSubMaxCSV;
  float fjMinCSV;
  float fjDoubleCSV;
  int fjgbb;
  int fjNbs;
  float fjGenPt;
  float fjGenSize;
  int fjIsMatched;
  float fjGenWPt;
  float fjGenWSize;
  int fjIsWMatched;
  int fjHighestPtGen;
  float fjHighestPtGenPt;
  int fjIsTight;
  int fjIsLoose;
  float fjRawPt;
  int fjNHF;
  float fjHTTMass;
  float fjHTTFRec;
  int fjIsClean;
  int fjNPartons;
  int fjNBPartons;
  int fjNCPartons;
  float fjPartonM;
  float fjPartonPt;
  float fjPartonEta;
  int fjGenNumB;
  int nHF;
  int nB;
  int nBGenJets;
  float genFatJetPt;
  int genFatJetNProngs;
  float fjsjPt[2];
  float fjsjEta[2];
  float fjsjPhi[2];
  float fjsjM[2];
  float fjsjCSV[2];
  float fjsjQGL[2];
  int nLoosePhoton;
  int nTightPhoton;
  int loosePho1IsTight;
  float loosePho1Pt;
  float loosePho1Eta;
  float loosePho1Phi;
  int loosePho1SelBit;
  int looseGenPho1PdgId;
  int nLooseLep;
  int nLooseElectron;
  int nLooseMuon;
  int nTightLep;
  int nTightElectron;
  int nTightMuon;
  float electronPt[4];
  float electronEta[4];
  float electronPhi[4];
  int electronSelBit[4];
  int electronPdgId[4];
  float electronSfLoose[4];
  float electronSfMedium[4];
  float electronSfTight[4];
  float electronSfMvaWP90[4];
  float electronSfMvaWP80[4];
  float electronSfUnc[4];
  float electronSfReco[4];
  float electronD0[4];
  float electronDZ[4];
  int electronNMissingHits[4];
  int electronTripleCharge[4];
  float electronCombIso[4];
  float muonPt[4];
  float muonEta[4];
  float muonPhi[4];
  int muonSelBit[4];
  int muonPdgId[4];
  float muonSfLoose[4];
  float muonSfMedium[4];
  float muonSfTight[4];
  float muonSfUnc[4];
  float muonSfReco[4];
  float muonD0[4];
  float muonDZ[4];
  int muonIsSoftMuon[4];
  float muonCombIso[4];
  float sf_zz;
  float sf_zzUnc;
  float sf_wz;
  float sf_vh;
  float sf_vhUp;
  float sf_vhDown;
  float genLep1Pt;
  float genLep1Eta;
  float genLep1Phi;
  int genLep1PdgId;
  float genLep2Pt;
  float genLep2Eta;
  float genLep2Phi;
  int genLep2PdgId;
  float genLep3Pt;
  float genLep3Eta;
  float genLep3Phi;
  int genLep3PdgId;
  float genLep4Pt;
  float genLep4Eta;
  float genLep4Phi;
  int genLep4PdgId;
  float genWPlusPt;
  float genWMinusPt;
  float genWPlusEta;
  float genWMinusEta;
  int looseGenLep1PdgId;
  int looseGenLep2PdgId;
  int looseGenLep3PdgId;
  int looseGenLep4PdgId;
  float diLepMass;
  int nTau;
  float mT[3];
  float hbbpt[3];
  float hbbeta[3];
  float hbbphi[3];
  float hbbm[3];
  float hbbm_reg[3];
  float hbbpt_reg[3];
  int nSoft2;
  int nSoft5;
  int nSoft10;
  float hbbCosThetaJJ[3];
  float hbbCosThetaCSJ1[3];
  float topMassLep1Met[3];
  float topWBosonCosThetaCS[3];
  float topWBosonPt;
  float topWBosonEta;
  float topWBosonPhi;
  float sumEtSoft1;
  float scaleUp;
  float scaleDown;
  float pdfUp;
  float pdfDown;
  float scale[6];
  float lheHT;
  int isGS;
};
#endif