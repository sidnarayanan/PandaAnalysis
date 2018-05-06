#ifndef MODULE
#define MODULE

#include "TString.h"
#include "vector"
#include "map"

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
        Container(T* ptr_) : ptr(ptr_) { }
        ~Container() { }
        T* ptr;
      };
      template <typename T>
      class ConstContainer : public BaseContainer {
      public:
        ConstContainer(const T* ptr_) : ptr(ptr_) { }
        ~ConstContainer() { }
        const T* ptr;
      };

    public:
      Registry() { }
      ~Registry() { for (auto& iter : _objs) { delete iter.second; } }
      template <typename T>
        void publish(TString name, T* ptr) { _objs[name] = new Container<T>(ptr); }
      template <typename T>
        void publishConst(TString name, const T* ptr) { _objs[name] = new ConstContainer<T>(ptr); }
      template <typename T>
        T* access(TString name) { 
          return dynamic_cast<Container<T>*>(_objs.at(name))->ptr; 
        }
      template <typename T>
        const T* accessConst(TString name) { 
          return dynamic_cast<ConstContainer<T>*>(_objs.at(name))->ptr; 
        }
      bool exists(TString name) { return _objs.find(name) != _objs.end(); }
    private:
      std::map<TString, BaseContainer*> _objs;
  };

  class BaseModule {
    public:
      BaseModule(TString name_): name(name_) { }
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

  class AnalysisMod : public BaseModule {
    public:
      AnalysisMod(TString name, 
                  panda::EventAnalysis& event_, 
                  Config& cfg_, 
                  Utils& utils_,
                  GeneralTree& gt_) : 
        BaseModule(name), 
        event(event_), 
        cfg(cfg_),
        utils(utils_),
        analysis(cfg.analysis),
        gt(gt_) { }
      virtual ~AnalysisMod() { for (auto* m : subMods) delete m; }
      
      // cascading calls to protected functions
      void initialize(Registry& registry);
      void readData(TString path);
      // execute DOES NOT cascade down child modules -
      // calling subMod execution is left up to the caller
      // to allow for more flexibility
      void execute();
      void reset();
      void terminate(); 

      virtual bool on() { return true; }

    protected:
      panda::EventAnalysis& event;
      Config& cfg;
      Utils& utils;
      const Analysis& analysis;
      GeneralTree& gt; 
      std::vector<AnalysisMod*> subMods; // memory management is done by parent

      // here, the module can publish and access data
      virtual void do_init(Registry& registry) { }
      virtual void do_readData(TString path) { }
      // this is where the actual execution is done
      virtual void do_execute() = 0;
      // reset objects between events, if needed
      virtual void do_reset() { }
      virtual void do_terminate() { }
  };
}

#endif
