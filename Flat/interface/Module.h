#ifndef MODULE
#define MODULE

#include "TString.h"
#include "vector"
#include "map"
#include "stdexcept"
#include "memory"

#include "Common.h"
#include "Config.h"

namespace pa {

  class Registry {
    private:
      class BaseContainer { // just for polymorphism
      public:
        BaseContainer() { }
        virtual ~BaseContainer() { }
      };
      template <typename T>
      class Container : public BaseContainer {
      public:
        Container(std::shared_ptr<T> ptr_) : ptr(ptr_) { }
        ~Container() { }
        const std::shared_ptr<T> ptr; // don't ask me why this is a copy and not a reference
                                      // when I tried using a reference the use_count went negative
                                      // ???
      };
      template <typename T>
        Container<T>* safe_cast(std::unique_ptr<BaseContainer>& base, TString name) {
          auto* cntr = dynamic_cast<Container<T>*>(base.get());
          if (cntr == nullptr) {
            PError("Registry::safe_cast", "Requesting object of wrong type: "+name+"!");
            throw std::runtime_error("");
          }
          return cntr;
        }

    public:
      Registry() { }
      ~Registry() { }
      template <typename T>
        void publish(TString name, std::shared_ptr<T>& ptr) {
          if (exists(name))
            PWarning("Registry::publish","UNDEFINED BEHAVIOR - multiple objects with name "+name+"!");
          _objs[name].reset(new Container<T>(ptr));
        }
      template <typename T>
        void publishConst(TString name, std::shared_ptr<T>& ptr) {
          if (exists(name))
            PWarning("Registry::publishConst","UNDEFINED BEHAVIOR - multiple objects with name "+name+"!");
          _const_objs[name].reset(new Container<T>(ptr));
        }
      template <typename T>
        std::shared_ptr<T> access(TString name, bool silentFail = false) {
          auto iter = _objs.find(name);
          if (iter == _objs.end()) {
            if (silentFail) {
              PWarning("Registry::access", "Could not access "+name+", returning (nil)!");
              return std::shared_ptr<T>(nullptr);
            } else {
              PError("Registry::access", "Could not access "+name+"!");
              throw std::runtime_error("");
            }
          }
          auto* cntr = safe_cast<T>(iter->second, name);
          return cntr->ptr;
        }
      template <typename T>
        std::shared_ptr<const T> accessConst(TString name, bool silentFail = false) {
          auto iter = _objs.find(name);
          if (iter == _objs.end()) {
            iter = _const_objs.find(name);
            if (iter == _const_objs.end()) {
              if (silentFail) {
                PWarning("Registry::accessConst", "Could not access "+name+", returning (nil)!");
                return std::shared_ptr<const T>(nullptr);
              } else {
                PError("Registry::accessConst", "Could not access "+name+"!");
                throw std::runtime_error("");
              }
            }
          }
          auto* cntr = safe_cast<T>(iter->second, name);
          return std::const_pointer_cast<const T>(cntr->ptr);
        }
      bool exists(TString name) { return !( _objs.find(name) == _objs.end()
                                            && _const_objs.find(name) == _const_objs.end() ); }
    private:
      // _objs can be accessed through access or accessConst; _const_objs only through accessConst
      std::map<TString, std::unique_ptr<BaseContainer>> _objs, _const_objs;
  };

  class BaseModule {
    public:
      BaseModule(TString name_): name(name_) {  }
      virtual ~BaseModule() { }

    protected:
      TString name;
  };

  class ConfigMod : public BaseModule {
    public:
      ConfigMod(const Analysis& a_, GeneralTree& gt, int DEBUG_);
      ~ConfigMod() { }

      void readData(TString path);
      const panda::utils::BranchList get_inputBranches() const { return bl; }

      Config cfg;
      Utils utils;

    protected:
      const Analysis& analysis;
      GeneralTree& gt;
      panda::utils::BranchList bl;

    private:
      void set_inputBranches();
      void set_outputBranches();
  };

  template <typename T>
  class BaseAnalysisMod : public BaseModule {
    public:
      BaseAnalysisMod(TString name,
                  panda::EventAnalysis& event_,
                  Config& cfg_,
                  Utils& utils_,
                  T& gt_,
                  int level_=0) :
        BaseModule(name),
        event(event_),
        cfg(cfg_),
        utils(utils_),
        analysis(cfg.analysis),
        gt(gt_),
        level(level_) { }
      virtual ~BaseAnalysisMod() {
        if (cfg.DEBUG > level + 2)
          PDebug("BaseAnalysisMod::~BaseAnalysisMod", name);
      }

      // cascading calls to protected functions
      void initialize(Registry& registry);
      void readData(TString path);
      // execute DOES NOT cascade down child modules -
      // calling subMod execution is left up to the caller
      // to allow for more flexibility
      void execute();
      void reset();
      void terminate();
      void print();

      virtual bool on() { return true; }

      template <typename MOD>
      MOD* addSubMod() {
        // add a sub module that takes a specific constructor signature
        auto* mod = new MOD(event, cfg, utils, gt, level + 1);
        subMods.emplace_back(mod);
        return mod;
      }

    protected:
      panda::EventAnalysis& event;
      Config& cfg;
      Utils& utils;
      const Analysis& analysis;
      T& gt;
      std::vector<std::unique_ptr<BaseAnalysisMod>> subMods; // memory management is done by parent
      int level;

      std::vector<TString> dump();
      // here, the module can publish and access data
      virtual void do_init(Registry& registry) { }
      virtual void do_readData(TString path) { }
      // this is where the actual execution is done
      virtual void do_execute() = 0;
      // reset objects between events, if needed
      virtual void do_reset() { }
      virtual void do_terminate() { }
  };
  typedef BaseAnalysisMod<GeneralTree> AnalysisMod; 
  typedef BaseAnalysisMod<HeavyResTree> HRMod; 

  // a completely empty mod
  class ContainerMod : public AnalysisMod {
    public:
      ContainerMod(TString name,
                   panda::EventAnalysis& event_,
                   Config& cfg_,
                   Utils& utils_,
                   GeneralTree& gt_,
                   int level_=0) :
        AnalysisMod(name, event_, cfg_, utils_, gt_, level_) { }
      ~ContainerMod() { }
      virtual bool on() { return true; }
    protected:
      virtual void do_execute() {
        for (auto& m : subMods)
          m->execute();
      }
      virtual void do_init(Registry& registry) {
        registry.access<std::vector<JESHandler>>("jesShifts");
      }
  };
}

#endif
