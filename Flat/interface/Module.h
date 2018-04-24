#include "TString.h"
#include "vector"
#include "map"


class Registry {
  public:
    Registry() { }
    ~Registry() { }
    template <typename T>
      void publish(TString name, T* ptr) { _objs[name] = static_cast<void*>(ptr); }
    template <typename T>
      T* access(TString name) { return static_cast<T*>(_objs.at(name)); }
  private:
    std::map<TString, void*> _objs;
};

class BaseModule {
  public:
    BaseModule(TString name_): name(name_) { }
    virtual ~BaseModule() { }

  private:
    TString name;
};

class ConfigMod : public BaseModule {
  public:
    ConfigMod(Analysis* a_, int DEBUG_) : 
      BaseModule("config") 
      cfg(a_, DEBUG_) { }
    ~ConfigMod() { }

    const Config& get_config() const { return cfg; }

  protected:
    Config cfg;
};

class AnalysisMod : public BaseModule {
  public:
    AnalysisMod(TString name, const panda::EventAnalysis& event_, const Config& cfg_, GeneralTree& gt_) : 
      BaseModule(name, cfg(cfg_)), event(event_), gt(gt_) { }
    ~AnalysisMod() { }
    
    // here, the module can publish and access data
    virtual initialize(Registry& registry) = 0;

    // this is where the actual execution is done
    virtual execute() = 0;

  protected:
    const Config& cfg;
    const panda::EventAnalysis& event;
    GeneralTree& gt; 
};
