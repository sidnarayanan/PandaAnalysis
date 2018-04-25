#include "Module.h"

namespace pa {
  class JetMod : public AnalysisMod {
  public: 
    JetMod(const panda::EventAnalysis& event_, 
           const Config& cfg_,                 
           const Utils& utils_,                
           GeneralTree& gt_) :                 
      AnalysisMod("jet", event_, cfg_, utils_, gt_) { }
    ~JetMod () { 
      delete ak4JERReader;
      for (auto& iter : ak4ScaleReader)
        delete iter.second;
      for (auto& iter : ak4UncReader)
        delete iter.second;
    }
    
  protected:
    void do_readData(TString path);
    void do_init(Registry& registry);
    void do_execute();  
  private:
    std::map<TString,FactorizedJetCorrector*> ak4ScaleReader; //!< calculate JES on the fly
    std::map<TString,JetCorrectionUncertainty*> ak4UncReader; //!< calculate JES unc on the fly
    JERReader *ak4JERReader{nullptr}; //!< fatjet jet energy resolution reader
    JetCorrectionUncertainty *uncReaderAK4  {nullptr};        
    FactorizedJetCorrector   *scaleReaderAK4{nullptr};        
    std::vector<JESHandler>* jesShifts{nullptr}; 

    void setupJES();
  };
}