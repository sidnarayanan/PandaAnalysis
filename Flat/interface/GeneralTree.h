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

    // public config
    bool monohiggs=false, vbf=false, fatjet=true, leptonic=false, photonic=false, hfCounting=false;
    bool btagWeights=false, useCMVA=false;
    std::map<ECFParams,float> fj1ECFNs;
    std::map<BTagParams,float> sf_btags;
    std::map<BTagParams,float> sf_alt_btags;
    std::map<TString,float> signal_weights;
    std::map<csvShift,float> sf_csvWeights;

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
  float pfmet[4];
  float pfmetphi;
  float pfmetnomu;
  float puppimet;
  float puppimetphi;
  float calomet;
  float calometphi;
  float pfcalobalance;
  float sumET;
  float trkmet;
  float trkmetphi;
  int whichRecoil;
  float puppiUWmag[4];
  float puppiUZmag[4];
  float puppiUAmag[4];
  float puppiUperp[4];
  float puppiUpara[4];
  float pfUWmag[4];
  float pfUZmag[4];
  float puppiUmag[4];
  float pfUAmag[4];
  float pfUperp[4];
  float pfUpara[4];
  float pfUmag[4];
  float puppiUWphi;
  float puppiUZphi;
  float puppiUAphi;
  float puppiUphi;
  float pfUWphi;
  float pfUZphi;
  float pfUAphi;
  float pfUphi;
  float dphipfmet[4];
  float dphipuppimet[4];
  float dphipuppiUW[4];
  float dphipuppiUZ[4];
  float dphipuppiUA[4];
  float dphipfUW[4];
  float dphipfUZ[4];
  float dphipfUA[4];
  float dphipuppiU[4];
  float dphipfU[4];
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
  int nJet[4];
  int nJot[4];
  int nIsoJet[4];
  float jetPt[2];
  float jetPEta[2];
  float jetPhi[2];
  float jetGenPt[2];
  float jetCSV[2];
  int jetFlav[2];
  int jetIsTight[2];
  int jetIsIso[2];
  float jotPt[4][NJET];
  float jotE[4][NJET];
  float jotEta[NJET];
  float jotPhi[NJET];
  float jotCSV[NJET];
  int jotVBFID[NJET];
  float jotCMVA[NJET];
  int jotIso[NJET];
  float jotQGL[NJET];
  float jotLep1Pt[NJET];
  float jotLep1PtRel[NJET];
  float jotLep1DeltaR[NJET];
  float jotTrk1Pt[NJET];
  float jotVtxPt[NJET];
  float jotVtxMass[NJET];
  float jotVtx3DVal[NJET];
  float jotVtx3DErr[NJET];
  int jotVtxNtrk[NJET];
  float jotEMF[NJET];
  float jotHF[NJET];
  int jotNLep[NJET];
  float jotGenPt[NJET];
  int jotFlav[NJET];
  int hbbjtidx[2];
  float jotRegFac[2];
  float barrelJet1Pt;
  float barrelJet1Eta;
  float barrelHT;
  float barrelHTMiss;
  float barrelJet12Pt;
  float jot12Mass[4];
  float jot12DEta[4];
  float jot12DPhi[4];
  int jetNBtags;
  int jetNMBtags;
  int isojetNBtags;
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
  int nHF;
  int nB;
  int nBGenJets;
  float genFatJetPt;
  int genFatJetNProngs;
  float fjsjPt[NSUBJET];
  float fjsjEta[NSUBJET];
  float fjsjPhi[NSUBJET];
  float fjsjM[NSUBJET];
  float fjsjCSV[NSUBJET];
  float fjsjQGL[NSUBJET];
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
  float electronPt[NLEP];
  float electronEta[NLEP];
  float electronPhi[NLEP];
  int electronSelBit[NLEP];
  int electronPdgId[NLEP];
  float electronSfLoose[NLEP];
  float electronSfMedium[NLEP];
  float electronSfTight[NLEP];
  float electronSfMvaWP90[NLEP];
  float electronSfMvaWP80[NLEP];
  float electronSfUnc[NLEP];
  float electronSfReco[NLEP];
  float electronD0[NLEP];
  float electronDZ[NLEP];
  int electronNMissingHits[NLEP];
  int electronTripleCharge[NLEP];
  float electronCombIso[NLEP];
  float muonPt[NLEP];
  float muonEta[NLEP];
  float muonPhi[NLEP];
  int muonSelBit[NLEP];
  int muonPdgId[NLEP];
  float muonSfLoose[NLEP];
  float muonSfMedium[NLEP];
  float muonSfTight[NLEP];
  float muonSfUnc[NLEP];
  float muonSfReco[NLEP];
  float muonD0[NLEP];
  float muonDZ[NLEP];
  int muonIsSoftMuon[NLEP];
  float muonCombIso[NLEP];
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
  int looseGenLep1PdgId;
  int looseGenLep2PdgId;
  int looseGenLep3PdgId;
  int looseGenLep4PdgId;
  float diLepMass;
  int nTau;
  float mT;
  float hbbpt[4];
  float hbbeta;
  float hbbphi;
  float hbbm[4];
  float hbbm_reg[4];
  int nSoft2;
  int nSoft5;
  int nSoft10;
  float hbbCosThetaJJ;
  float hbbCosThetaCSJ1;
  float topMassLep1Met[4];
  float topWBosonCosThetaCS;
  float topWBosonPt;
  float topWBosonEta;
  float topWBosonPhi;
  float sumEtSoft1;
  float scaleUp;
  float scaleDown;
  float pdfUp;
  float pdfDown;
  float scale[6];
  int isGS;
};
#endif