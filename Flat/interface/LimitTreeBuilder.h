
#ifndef LimitTreeBuilder_h
#define LimitTreeBuilder_h

#include <TROOT.h>
#include <TTree.h>
#include <TTreeFormula.h>
#include <TFile.h>
#include <TSelector.h>
#include <TString.h>
#include <map>
#include <vector>
#include "PandaCore/Tools/interface/Common.h"

// used to give distinct names to trees
extern int treeCounter;

class xformula {
	public:
		xformula(TString n_, TString f_) { name=n_; formula=f_; val = new float(0); }
		~xformula() { delete val;	}
		TString name;
		TString formula;
		float *val;
};


// used to map input variables
class VariableMap {
public:
	VariableMap() {
	}
	~VariableMap() {}
	void AddVar(TString n,TString v) { vmap[n] = v; }
	void AddFormula(TString n,TString f) { fmap[n] = f; }
	const std::map<TString,TString> GetVars() const { return vmap; }
	const std::map<TString,TString> GetFormulae() const { return fmap; }
private:
	std::map<TString,TString> vmap, fmap;
};

// defines a process for a given region and choice of systematic shifts
class Process {
public:
	Process(TString n, TTree *in, VariableMap *vPtr, TString sel, TString w);
	~Process();
	void Run();
	TTree *GetTree() { return limitTree; }
	TTree *GetInput() { return inputTree; }
	TString name;
	TString syst="";
private:
	TTree *limitTree=0;
	TTree *inputTree=0;
	TString selection;
	TString weight;
	std::vector<xformula*> vars, formulae;
};

// defines a region, which is essentially just a list of Processes
class Region {
public:
	Region(TString n) { name = n; }
	~Region() {}
	void AddProcess(Process *p) { ps.push_back(p); }
	std::vector<Process*> GetProcesses() { return ps; }
	void Run() { 
		for (auto p : ps) { 
			PInfo("LimitTreeBuilder::Region::Run",TString::Format("%s",name.Data())); 
			p->Run(); 
		} 
	}
	TString name;
private:
	std::vector<Process*> ps;
};

class LimitTreeBuilder {
public :

	LimitTreeBuilder() {}
	~LimitTreeBuilder() {}
	void SetOutFile(TString f) { fOut = new TFile(f,"RECREATE"); }
	void AddRegion(Region *r) { regions.push_back(r); }
	void cd() { fOut->cd(); }
	void Run() { for (auto r : regions) { fOut->cd(); r->Run(); } }
	void Output();
private:
	std::vector<Region*> regions;
	TFile *fOut=0;
};

#endif

