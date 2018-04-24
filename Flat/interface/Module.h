#include "TString.h"
#include "vector"
#include "map"

namespace pa {

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

      void setDataDir(TString path);

      const Config& get_config() const { return cfg; }
      const Utils& get_utils() const { return utils; }

    protected:
      const Analysis& analysis;
      Config cfg;
      Utils utils;
  };

  class AnalysisMod : public BaseModule {
    public:
      AnalysisMod(TString name, 
                  const panda::EventAnalysis& event_, 
                  const Config& cfg_, 
                  const Utils& utils_,
                  GeneralTree& gt_) : 
        BaseModule(name), 
        event(event_), 
        cfg(cfg_),
        utils(utils_),
        gt(gt_) { }
      virtual ~AnalysisMod() { }
      
      // here, the module can publish and access data
      virtual initialize(Registry& registry) = 0;

      // this is where the actual execution is done
      virtual execute() = 0;

      // reset objects between events, if needed
      virtual reset() { }

    protected:
      const panda::EventAnalysis& event;
      const Config& cfg;
      const Utils& utils;
      GeneralTree& gt; 
  };
}

