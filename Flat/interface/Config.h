#include "TString.h"
#include "vector"
#include "map"
// fastjet
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/MeasureDefinition.hh"

// btag
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
//#include "BTagCalibrationStandalone.h"

// JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// Utils
#include "PandaAnalysis/Utilities/interface/RoccoR.h"
#include "PandaAnalysis/Utilities/interface/CSVHelper.h"
#include "PandaAnalysis/Utilities/interface/EnergyCorrelations.h"

namespace pa {
  enum CorrectionType { //!< enum listing relevant corrections applied to MC
      cNPV=0,       //!< npv weight
      cPU,          //!< true pu weight
      cPUUp,        //!< true pu weight
      cPUDown,      //!< true pu weight
      cEleVeto,     //!< monojet SF, Veto ID for e
      cEleLoose,    //!< monojet SF, Tight ID for e
      cEleMedium,   //!< monojet SF, Tight ID for e
      cEleTight,    //!< monojet SF, Tight ID for e
      cEleMvaWP90,
      cEleMvaWP80,
      cEleReco,     //!< monojet SF, tracking for e
      cWmHEwkCorr,     //!< W(l-V)H Ewk Corr weight  
      cWmHEwkCorrUp,   //!< W(l-V)H Ewk Corr weight Up  
      cWmHEwkCorrDown, //!< W(l-V)H Ewk Corr weight Down  
      cWpHEwkCorr,     //!< W(l+v)H Ewk Corr weight  
      cWpHEwkCorrUp,   //!< W(l+v)H Ewk Corr weight Up  
      cWpHEwkCorrDown, //!< W(l+v)H Ewk Corr weight Down  
      cZnnHEwkCorr,     //!< Z(vv)H Ewk Corr weight  
      cZnnHEwkCorrUp,   //!< Z(vv)H Ewk Corr weight Up  
      cZnnHEwkCorrDown, //!< Z(vv)H Ewk Corr weight Down  
      cZllHEwkCorr,     //!< Z(ll)H Ewk Corr weight  
      cZllHEwkCorrUp,   //!< Z(ll)H Ewk Corr weight Up  
      cZllHEwkCorrDown, //!< Z(ll)H Ewk Corr weight Down  
      cWZEwkCorr,
      cqqZZQcdCorr,
      cMuLooseID,   //!< MUO POG SF, Loose ID for mu 
      cMuMediumID,  //!< MUO POG SF, Tight ID for mu 
      cMuTightID,   //!< MUO POG SF, Tight ID for mu 
      cMuLooseIso,  //!< MUO POG SF, Loose Iso for mu 
      cMuMediumIso, //!< MUO POG SF, Loose Iso for mu 
      cMuTightIso,  //!< MUO POG SF, Tight Iso for mu 
      cMuReco,      //!< MUO POG SF, tracking for mu
      cPho,         //!< EGM POG SF, contains ID for gamma
      cTrigMET,     //!< MET trigger eff        
      cTrigMETZmm,  //!< Zmumu MET trigger eff
      cTrigEle,     //!< Ele trigger eff        
      cTrigMu,      //!< Mu trigger eff        
      cTrigPho,     //!< Pho trigger eff        
      cZNLO,        //!< NLO weights for QCD Z,W,A,A+2j
      cWNLO,
      cANLO,
      cANLO2j,
      cZEWK,        //!< EWK weights for QCD Z,W,A,A+2j
      cWEWK,
      cAEWK,
      cVBF_ZNLO,    //!< NLO weights for QCD Z,W in VBF phase space
      cVBF_ZllNLO,  
      cVBF_WNLO,
      cVBFTight_ZNLO,    //!< NLO weights for QCD Z,W in tight VBF phase space
      cVBFTight_ZllNLO,  
      cVBFTight_WNLO,
      cVBF_EWKZ,    //!< k-factors for EWK Z,W in VBF phase space
      cVBF_EWKW,
      cVBF_TrigMET, //!< MET trigger eff as a f'n of mjj/met 
      cVBF_TrigMETZmm,
      cBadECALJets,  //!< bad ECAL clusters to filter jets
      cJetLoosePUID,
      cN
  };

  enum LepSelectionBit {
   kLoose     =(1<<0),
   kFake      =(1<<1),
   kMedium    =(1<<2),
   kTight     =(1<<3),
   kDxyz      =(1<<4),
   kEleMvaWP90=(1<<5),
   kEleMvaWP80=(1<<6)
  };

  enum PhoSelectionBit {
   pMedium    =(1<<0),
   pTight     =(1<<1),
   pHighPt    =(1<<2),
   pCsafeVeto =(1<<3),
   pPixelVeto =(1<<4),
   pTrkVeto   =(1<<5)
  };

  enum BTagType {
      bJetL=0,
      bSubJetL,
      bJetM,
      bSubJetM,
      bN
  };

  enum TriggerBits {
      kMETTrig       = 0,
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

  struct GenJetInfo {
    public:
      float pt=-1, eta=-1, phi=-1, m=-1;
      float msd=-1;
      float tau3=-1, tau2=-1, tau1=-1;
      float tau3sd=-1, tau2sd=-1, tau1sd=-1;
      int nprongs=-1;
      float partonpt=-1, partonm=-1;
      std::vector<std::vector<float>> particles;
      std::vector<std::vector<std::vector<float>>> ecfs; // uh
      void reset() {
        pt=-1; eta=-1; phi=-1; m=-1;
        msd=-1;
        tau3=-1; tau2=-1; tau1=-1;
        tau3sd=-1; tau2sd=-1; tau1sd=-1;
        nprongs=-1;
        partonpt=-1; partonm=-1;
        for (auto& v : particles) {
          std::fill(v.begin(), v.end(), 0);
        }
        for (auto& v : ecfs) {
          for (auto& vv : v) {
            std::fill(vv.begin(), vv.end(), -1);
          }
        }
      }
  };

  struct JetHistory {
    int user_idx;
    int child_idx;
  };

  struct JetWrapper {
      float pt; 
      int flavor{0};
      float genpt{0};
      bool iso{false};
      const panda::Jet* base;

      JetWrapper(float pt_, const panda::Jet& j): pt(pt_), base(&j) { }

      const panda::Jet& get_base() const { return *base; }
      float scale() const { return pt / base->pt(); }
      TLorentzVector p4() const { 
          TLorentzVector v; 
          v.SetPtEtaPhiM(pt, base->eta(), base->phi(), base->m()); 
          return v;
      }
  };

  struct JESHandler {
      std::vector<JetWrapper> all; // all jets 
      std::vector<const JetWrapper*> cleaned;
      std::vector<const JetWrapper*> iso;     // cleaned that do not overlap with fj
      std::vector<const JetWrapper*> central; // cleaned that are central
      std::vector<const JetWrapper*> bcand;   // all that are b candidates

      TLorentzVector vpfMET, vpuppiMET;
      TVector2 vpfMETNoMu;
      TLorentzVector vpfUW, vpfUZ, vpfUA;
      TLorentzVector vpuppiUW, vpuppiUZ, vpuppiUA;

      void clear() { 
              all.clear(); cleaned.clear(); iso.clear();
              central.clear(); bcand.clear(); 
              vpfMETNoMu.SetMagPhi(0,0);
              for (TLorentzVector* v_ : {&vpfMET, &vpuppiMET, 
                                         &vpfUW, &vpfUZ, &vpfUA,
                                         &vpuppiUW, &vpuppiUZ, &vpuppiUA})
                  v_->SetPtEtaPhiM(0,0,0,0);
          }
      void reserve(int N) { all.reserve(N); cleaned.reserve(N); iso.reserve(N);
                            central.reserve(N); bcand.reserve(N); }
  };

  class Config {
    public:
      Config(const Analysis& a_, int DEBUG_ = 0) : 
        DEBUG(DEBUG_),
        analysis(a_),
        tr("PandaAnalyzer::Run", DEBUG+1)
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

      bool isData{false};              // to do gen matching, etc

      int NPFPROPS{9}, NSVPROPS{13}, NMAXPF{100}, NMAXSV{10}, NGENPROPS{8};
      float FATJETMATCHDR2{2.25};


      // configuration read from output tree
      std::vector<int> ibetas;
      std::vector<int> Ns; 
      std::vector<int> orders;
  };

  class Utils {
    public:
      Utils() { }
      ~Utils(); 

      double getCorr(CorrectionType ct, double x, double y=0);
      double getError(CorrectionType ct, double x, double y=0);
      void openCorr(CorrectionType, TString, TString, int);

      std::map<int, int> pdgToQ; 
      
      // fastjet reclustering
      fastjet::JetDefinition       *jetDef                {nullptr};
      fastjet::JetDefinition       *jetDefKt              {nullptr};
      fastjet::contrib::SoftDrop   *softDrop              {nullptr};
      fastjet::contrib::Njettiness *tauN                  {nullptr};
      fastjet::AreaDefinition      *areaDef               {nullptr};
      fastjet::GhostedAreaSpec     *activeArea            {nullptr};
      fastjet::JetDefinition       *jetDefGen             {nullptr};
      fastjet::JetDefinition       *softTrackJetDefinition{nullptr};

      // CMSSW-provided utilities
      BTagCalibration *btagCalib   {nullptr};
      BTagCalibration *sj_btagCalib{nullptr};
      std::vector<BTagCalibrationReader*> btagReaders = std::vector<BTagCalibrationReader*>(bN,0); 
        //!< maps BTagType to a reader 

      Binner btagpt  = Binner({});
      Binner btageta = Binner({});
      std::vector<std::vector<double>> lfeff, ceff, beff;
      TMVA::Reader *bjetregReader{nullptr}; 
      float *bjetreg_vars{nullptr};

      std::map<TString,JetCorrectionUncertainty*> ak8UncReader; //!< calculate JES unc on the fly
      JERReader *ak8JERReader{nullptr}; //!< fatjet jet energy resolution reader

      JetCorrectionUncertainty *uncReader     {nullptr};           

      EraHandler eras = EraHandler(2016); //!< determining data-taking era, to be used for era-dependent JEC
      ParticleGridder   *grid   {nullptr};
      pa::ECFCalculator *ecfcalc{nullptr};

      auto   fCorrs = std::vector<TFile*>    (cN,nullptr);
      auto  f1Corrs = std::vector<TF1Corr*>  (cN,nullptr);
      auto  h1Corrs = std::vector<TH1Corr*>  (cN,nullptr);
      auto  h2Corrs = std::vector<TH2Corr*>  (cN,nullptr);

      TFile     *fMSDcorr             {nullptr};
      TF1       *puppisd_corrGEN      {nullptr};
      TF1       *puppisd_corrRECO_cen {nullptr};
      TF1       *puppisd_corrRECO_for {nullptr};
      RoccoR    *rochesterCorrection  {nullptr};
      CSVHelper *csvReweighter        {nullptr},   *cmvaReweighter  {nullptr};

  };

}

