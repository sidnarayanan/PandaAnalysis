#ifndef OPERATOR
#define OPERATOR

#include "TString.h"
#include "vector"
#include "map"
#include "stdexcept"
#include "memory"

#include "PandaAnalysis/Flat/interface/Common.h"
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
        // default copy and assignment are fine 
        const std::shared_ptr<T> ptr; // don't ask me why this is a copy and not a reference
                                      // when I tried using a reference the use_count went negative
                                      // ???
      };

      template <typename T>
      Container<T>* safe_cast(std::shared_ptr<BaseContainer>& base, TString name) {
        auto* cntr = dynamic_cast<Container<T>*>(base.get());
        if (cntr == nullptr) {
          logger.error("Registry::safe_cast", "Requesting object of wrong type: "+name+"!");
          throw std::runtime_error("");
        }
        return cntr;
      }

    public:
      Registry() { }
      ~Registry() { }
      Registry(const Registry& other) : _objs(other._objs), _const_objs(other._const_objs) { }
      Registry& operator=(const Registry& other) { 
        if (&other == this)
          return *this; 
        _objs = other._objs;
        _const_objs = other._const_objs;
        return *this; 
      }

      template <typename T>
      void publish(TString name, std::shared_ptr<T>& ptr) {
        if (exists(name))
          logger.warning("Registry::publish","Undefined behavior - multiple objects with name "+name+"!");
        _objs[name].reset(new Container<T>(ptr));
      }

      template <typename T>
      void publishConst(TString name, std::shared_ptr<T>& ptr) {
        if (exists(name))
          logger.warning("Registry::publishConst","Undefined behavior - multiple objects with name "+name+"!");
        _const_objs[name].reset(new Container<T>(ptr));
      }

      template <typename T>
      std::shared_ptr<T> access(TString name, bool silentFail = false) {
        auto iter = _objs.find(name);
        if (iter == _objs.end()) {
          if (silentFail) {
            logger.warning("Registry::access", "Could not access "+name+", returning (nil)!");
            return std::shared_ptr<T>(nullptr);
          } else {
            logger.error("Registry::access", "Could not access "+name+"!");
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
              logger.warning("Registry::accessConst", "Could not access "+name+", returning (nil)!");
              return std::shared_ptr<const T>(nullptr);
            } else {
              logger.error("Registry::accessConst", "Could not access "+name+"!");
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
      std::map<TString, std::shared_ptr<BaseContainer>> _objs, _const_objs;
  };
  
  template <typename T>
  class BaseOperator {
    public:
      BaseOperator(TString name_, T& gt_): name(name_), gt(gt_) {  }
      virtual ~BaseOperator() { }
      BaseOperator(const BaseOperator&) = delete;
      BaseOperator& operator=(const BaseOperator&) = delete; 

    protected:
      TString name;
      T& gt; 
  };

  class ConfigOp : public BaseOperator<GeneralTree> {
    public:
      ConfigOp(Analysis& a_, GeneralTree& gt, int DEBUG_);
      ~ConfigOp() { }

      void readData(TString path);
      const panda::utils::BranchList get_inputBranches() const { return bl; }

      Config cfg;
      Utils utils;

    protected:
      Analysis& analysis;
      panda::utils::BranchList bl;

    private:
      void set_inputBranches();
      void set_outputBranches();
  };

  template <typename T>
  class BaseAnalysisOp : public BaseOperator<T> {
    public:
      BaseAnalysisOp(TString name_,
                      panda::EventAnalysis& event_,
                      Config& cfg_,
                      Utils& utils_,
                      T& gt_,
                      int level_=0) :
        BaseOperator<T>(name_, gt_),
        event(event_),
        cfg(cfg_),
        utils(utils_),
        analysis(cfg.analysis),
        level(level_) { }
      virtual ~BaseAnalysisOp() {
        if (cfg.DEBUG > level + 2)
          logger.debug("BaseAnalysisOp::~BaseAnalysisOp", this->name);
      }

      // cascading calls to protected functions
      void initialize(Registry& registry);
      void readData(TString path);
      // execute DOES NOT cascade down child operators -
      // calling subOp execution is left up to the caller
      // to allow for more flexibility
      void execute();
      void reset();
      void terminate();
      void print();

      virtual bool on() { return true; }

      template <typename OP>
      OP* addSubOp() {
        // add a sub operator that takes a specific constructor signature
        auto* op = new OP(event, cfg, utils, this->gt, level + 1);
        subOps.emplace_back(op);
        return op;
      }

    protected:
      panda::EventAnalysis& event;
      Config& cfg;
      Utils& utils;
      const Analysis& analysis;
      std::vector<std::unique_ptr<BaseAnalysisOp>> subOps; // memory management is done by parent
      int level;

      std::vector<TString> dump();
      // here, the operator can publish and access data
      virtual void do_init(Registry& registry) { }
      virtual void do_readData(TString path) { }
      // this is where the actual execution is done
      virtual void do_execute() = 0;
      // reset objects between events, if needed
      virtual void do_reset() { }
      virtual void do_terminate() { }
  };
  typedef BaseAnalysisOp<GeneralTree> AnalysisOp; 
  typedef BaseAnalysisOp<HeavyResTree> HROp; 

  // a completely empty op
  class ContainerOp : public AnalysisOp {
    public:
      ContainerOp(TString name,
                   panda::EventAnalysis& event_,
                   Config& cfg_,
                   Utils& utils_,
                   GeneralTree& gt_,
                   int level_=0) :
        AnalysisOp(name, event_, cfg_, utils_, gt_, level_) { }
      ~ContainerOp() { }
      virtual bool on() { return true; }
    protected:
      virtual void do_execute() {
        for (auto& m : subOps)
          m->execute();
      }
      virtual void do_init(Registry& registry) {
        registry.access<std::vector<JESHandler>>("jesShifts");
      }
  };
}

#endif
