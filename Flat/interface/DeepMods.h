#ifndef DEEPMODS
#define DEEPMODS

#include "Module.h"
#include "set"
#include "PandaAnalysis/Utilities/interface/EnergyCorrelations.h"
#include "fastjet/contrib/Njettiness.hh"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

namespace pa {

  /*
  class BRegBDTMod : public Analysis {
  public:
    BRegBDTMod(panda::EventAnalysis& event_,
                Config& cfg_,
                Utils& utils_,
                GeneralTree& gt_) :
      AnalysisMod("bregbdt", event_, cfg_, utils_, gt_) { 
      }

    virtual bool on() { return !analysis.genOnly && analysis.hbb && analysis.bjetRegression; }
  protected:
    void do_readData(TString dirPath) {
      TString modelfile = dirPath+"/trainings/breg_training_2017.pb";
      downloadData("http://t3serv001.mit.edu/~snarayan/pandadata/trainings/breg_training_2017.pb",
                   modelfile);
      build(modelfile);
    }
    virtual void do_init(Registry& registry) {
      TFInferMod::do_init(registry);
      currentJet = registry.access<JetWrapper*>("higgsDaughterJet");
    }
    void do_execute(); 
  private:
    JetWrapper **currentJet{nullptr};
  };

    */

  class TFInferMod : public AnalysisMod {
  public:
    TFInferMod(TString name_,
               panda::EventAnalysis& event_,
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_) :
      AnalysisMod(name_, event_, cfg_, utils_,  gt_),
      outputNames(1) { }
    ~TFInferMod() { }
    
  protected:
    virtual void do_init(Registry& registry) {
      // just in case you want to set these in other modules 
      registry.publishConst(name+"_outputs", &outputs);
      registry.publish(name+"_inputs", &inputs);
    }
    virtual void do_readData(TString dirPath) = 0; // must call build(weightpath) 
    virtual void do_execute() = 0; // must call eval() 
    virtual void do_reset() { 
      fill(inputs.begin(), inputs.end(), -99); 
      fill(outputs.begin(), outputs.end(), -99); 
    }

    virtual void build(TString weightpath) final;
    virtual void eval() final;

    tensorflow::GraphDef* graph{nullptr}; // need to figure out if these are owned by tf or not
    tensorflow::Session* sess{nullptr};
    TString inputName, outputName;
    int n_inputs{0}, n_outputs{0}; 
    std::vector<float> inputs, outputs;
  private:
    std::vector<std::string> outputNames; 
  };

  class BRegDeepMod : public TFInferMod {
  public:
    BRegDeepMod(panda::EventAnalysis& event_,
                Config& cfg_,
                Utils& utils_,
                GeneralTree& gt_) :
      TFInferMod("bregdeep", event_, cfg_, utils_, gt_) { 
        n_inputs = 43;
        n_outputs = 3; 
        inputName = "ffwd_inp";
        outputName = "ffwd_out/BiasAdd";
      }

    virtual bool on() { return !analysis.genOnly && analysis.hbb && analysis.bjetRegression; }
  protected:
    void do_readData(TString dirPath) {
      TString modelfile = dirPath+"/trainings/breg_training_2017.pb";
      downloadData("http://t3serv001.mit.edu/~snarayan/pandadata/trainings/breg_training_2017.pb",
                   modelfile);
      build(modelfile);
    }
    virtual void do_init(Registry& registry) {
      TFInferMod::do_init(registry);
      currentJet = registry.access<JetWrapper*>("higgsDaughterJet");
    }
    void do_execute(); 
  private:
    JetWrapper **currentJet{nullptr};
  };

  template <typename GENP>
  class DeepGenMod : public AnalysisMod {
  public: 
    DeepGenMod(panda::EventAnalysis& event_, 
               Config& cfg_,
               Utils& utils_,
               GeneralTree& gt_);
    virtual ~DeepGenMod () { 
      delete jetDef; delete tauN;
      delete ecfcalc; delete grid;
    }

    virtual bool on() { return !analysis.isData && analysis.deepGen; }
    
  protected:
    void do_init(Registry& registry) {
      genP = registry.accessConst<std::vector<panda::Particle*>>("genP");
      fOut = registry.access<TFile>("fOut");
      incrementAux(false); 
    }
    void do_execute(); 
    void do_reset() { 
      genJetInfo.reset();
      if (grid != nullptr)
        grid->clear();
    }
    void do_terminate() {
      incrementAux(true);
    }
    void countGenPartons(std::unordered_set<const GENP*>&);
    void incrementAux(bool close = false);
  private:
    const std::vector<panda::Particle*> *genP{nullptr};
    fastjet::JetDefinition       *jetDef                {nullptr};
    fastjet::contrib::Njettiness *tauN                  {nullptr};
    ECFCalculator *ecfcalc{nullptr};
    ParticleGridder *grid{nullptr};
    GenJetInfo genJetInfo;
    TFile *fAux{nullptr};
    TTree *tAux{nullptr};
    TFile *fOut{nullptr};
    int auxCounter{0}; 
  };

  typedef DeepGenMod<panda::GenParticle> DeepPGenMod;
  typedef DeepGenMod<panda::UnpackedGenParticle> DeepUGenMod;
}


#endif
