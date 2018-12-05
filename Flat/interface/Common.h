#ifndef PACOMMON
#define PACOMMON

// STL
#include "vector"
#include <unordered_set>
#include "map"
#include <string>
#include <cmath>
#include "memory"

// ROOT
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>

#include "GeneralTree.h"
#include "HeavyResTree.h"
#include "L1Tree.h" 
#include "PandaTree/Objects/interface/EventAnalysis.h"
#include "PandaCore/Tools/interface/Common.h"
#include "PandaCore/Tools/interface/DataTools.h"
#include "PandaCore/Tools/interface/JERReader.h"

// macros
#define JESLOOP for (int shift = 0; shift != cfg.maxshiftJES; ++shift)
#define METLOOP for (int shift = 0; shift != 3; ++shift)

namespace pa {
  inline int jes2i(shiftjes i) { return static_cast<int>(i); }
  inline shiftjes i2jes(int i) { return static_cast<shiftjes>(i); }

  enum ProcessType { 
      kNoProcess,
      kZ,
      kW,
      kA,
      kZEWK,
      kWEWK,
      kTT,
      kTop, // used for non-ttbar top
      kV, // used for non V+jets W or Z
      kH,
      kSignal,
      kVV
  };
  
  class Analysis {
  public:
    Analysis(TString name_ = "") { name = name_; }
    ~Analysis() {}
    TString name;
    TString inpath="";
    TString outpath="";
    TString datapath="";
    ProcessType processType=kNoProcess;
    int year=2016;
    bool isData;
    
    bool ak = false; 
    bool ak8 = false;
    bool applyJER = false; 
    bool bjetRegTraining = false;
    bool bjetBDTReg = false;
    bool bjetDeepReg = false; 
    bool btagSFs = true;
    bool btagWeights = false;
    bool complicatedLeptons = false;
    bool complicatedPhotons = false;
    bool deep = false;
    bool deepAntiKtSort = false;
    bool deepExC = false;
    bool deepGen = false;
    bool deepGenGrid = false;
    bool deepKtSort = false;
    bool deepSVs = false;
    bool deepTracks = false;
    bool fatjet = true;
    bool firstGen = true;
    bool genOnly = false;
    bool hbb = false;
    bool hfCounting = false;
    bool jetFlavorPartons = true;
    bool jetFlavorJets = false;
    bool mcTriggers = false;
    bool monoh = false;
    bool puppiJets = true;
    bool recalcECF = false; 
    bool recluster = false;
    bool reclusterFJ = false;
    bool recoil = true;
    bool rerunJER = false;
    bool rerunJES = false;
    bool useCMVA = false;
    bool useDeepCSV = false;
    bool unpackedGen = false; 
    bool varyJES = false;
    bool varyJESTotal = false;
    bool vbf = false;
    bool vbfhbb = false; 
    bool vqqhbb = false; 
    bool zllhbb = false; 
  };

  enum CorrectionType { //!< enum listing relevant corrections applied to MC
    cNPV=0,     //!< npv weight
    cPU,      //!< true pu weight
    cPUUp,    //!< true pu weight
    cPUDown,    //!< true pu weight
    cEleVeto,   //!< monojet SF, Veto ID for e
    cEleLoose,  //!< monojet SF, Tight ID for e
    cEleMedium,   //!< monojet SF, Tight ID for e
    cEleTight,  //!< monojet SF, Tight ID for e
    cEleMvaWP90,
    cEleMvaWP80,
    cEleReco,   //!< monojet SF, tracking for e
    cWZEwkCorr,
    cqqZZQcdCorr,
    cMuLooseID,   //!< MUO POG SF, Loose ID for mu 
    cMuMediumID,  //!< MUO POG SF, Tight ID for mu 
    cMuTightID,   //!< MUO POG SF, Tight ID for mu 
    cMuLooseIso,  //!< MUO POG SF, Loose Iso for mu 
    cMuMediumIso, //!< MUO POG SF, Loose Iso for mu 
    cMuTightIso,  //!< MUO POG SF, Tight Iso for mu 
    cMuReco,    //!< MUO POG SF, tracking for mu
    cPho,     //!< EGM POG SF, contains ID for gamma
    cPhoFake,   //!< jet-faking-photon rate 
    cTrigMET,   //!< MET trigger eff    
    cTrigMETZmm,  //!< Zmumu MET trigger eff
    cTrigEle,   //!< Ele trigger eff    
    cTrigMu,    //!< Mu trigger eff    
    cTrigPho,   //!< Pho trigger eff    
    cTrigDoubleEleLeg1, //!< Double Ele trigger eff    
    cTrigDoubleEleLeg2, //!< Double Ele trigger eff    
    cTrigDoubleMuLeg1, //!< Double Mu trigger eff
    cTrigDoubleMuLeg2, //!< Double Mu trigger eff
    cZNLO,    //!< NLO weights for QCD Z,W,A,A+2j
    cWNLO,
    cANLO,
    cANLO2j,
    cZEWK,    //!< EWK weights for QCD Z,W,A,A+2j
    cWEWK,
    cAEWK,
    cVBF_ZNLO,  //!< NLO weights for QCD Z,W in VBF phase space
    cVBF_ZllNLO,  
    cVBF_WNLO,
    cVBFTight_ZNLO,  //!< NLO weights for QCD Z,W in tight VBF phase space
    cVBFTight_ZllNLO,  
    cVBFTight_WNLO,
    cVBF_EWKZ,  //!< k-factors for EWK Z,W in VBF phase space
    cVBF_EWKW,
    cVBF_TrigMET, //!< MET trigger eff as a f'n of mjj/met 
    cVBF_TrigMETZmm,
    cBadECALJets,  //!< bad ECAL clusters to filter jets
    cJetLoosePUID,
    cCSVBL, //!< CSV loose WP efficiencies 
    cCSVCL,
    cCSVLL,
    cWmHEwkCorr,   //!< W(l-V)H Ewk Corr weight  
    cWmHEwkCorrUp,   //!< W(l-V)H Ewk Corr weight Up  
    cWmHEwkCorrDown, //!< W(l-V)H Ewk Corr weight Down  
    cWpHEwkCorr,   //!< W(l+v)H Ewk Corr weight  
    cWpHEwkCorrUp,   //!< W(l+v)H Ewk Corr weight Up  
    cWpHEwkCorrDown, //!< W(l+v)H Ewk Corr weight Down  
    cZnnHEwkCorr,   //!< Z(vv)H Ewk Corr weight  
    cZnnHEwkCorrUp,   //!< Z(vv)H Ewk Corr weight Up  
    cZnnHEwkCorrDown, //!< Z(vv)H Ewk Corr weight Down  
    cZllHEwkCorr,   //!< Z(ll)H Ewk Corr weight  
    cZllHEwkCorrUp,   //!< Z(ll)H Ewk Corr weight Up  
    cZllHEwkCorrDown, //!< Z(ll)H Ewk Corr weight Down  
    cL1PreFiring, //!< PreFiring weights 
    cL1PhotonPreFiring, //!< PhotonPreFiring weights 
    cN
  };

  enum LepSelectionBit {
   kLoose   =(1<<0),
   kFake    =(1<<1),
   kMedium  =(1<<2),
   kTight   =(1<<3),
   kDxyz    =(1<<4),
   kEleMvaWP90=(1<<5),
   kEleMvaWP80=(1<<6),
   kMvaMedium    =(1<<7),
   kMvaTight     =(1<<8),
   kMiniIsoMedium=(1<<9),
   kMiniIsoTight =(1<<10)
  };

  enum PhoSelectionBit {
   pMedium  =(1<<0),
   pTight   =(1<<1),
   pHighPt  =(1<<2),
   pCsafeVeto =(1<<3),
   pPixelVeto =(1<<4),
   pTrkVeto   =(1<<5)
  };

  enum TriggerBits {
    kMETTrig     = 0,
    kSingleEleTrig,
    kSinglePhoTrig,
    kSingleMuTrig,
    kDoubleMuTrig,
    kDoubleEleTrig,
    kEMuTrig,
    kJetHTTrig,
    kMuFakeTrig,
    kEleFakeTrig,
    kNTrig
  };


}

#endif
