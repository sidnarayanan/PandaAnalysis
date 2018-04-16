#include <TTreeFormula.h>
#include "../interface/LimitTreeBuilder.h"
#include "PandaCore/Tools/interface/TreeTools.h"

Process::Process(TString n, TTree *in, VariableMap *vPtr, TString sel, TString w) {
	name = n;
	inputTree = in;
	selection = sel;
	weight = w;

	VariableMap v = *vPtr;
	limitTree = new TTree(TString::Format("tree%i",treeCounter++),TString::Format("limit tree for %s%s",n.Data(),syst.Data()));
	std::map<TString,TString> vars_ = vPtr->GetVars();
	for (auto it=vars_.begin(); it!=vars_.end(); ++it) {
		xformula *x = new xformula(it->first,it->second);
		vars.push_back(x);
		limitTree->Branch(x->name,x->val,x->name+"/F");
	}
	vars_ = vPtr->GetFormulae();
	for (auto it=vars_.begin(); it!=vars_.end(); ++it) {
		xformula *x = new xformula(it->first,it->second);
		formulae.push_back(x);
		limitTree->Branch(x->name,x->val,x->name+"/F");
	}
}

Process::~Process() {
	for (auto *x : vars)
		delete x;
	for (auto *x : formulae) 
		delete x;
	delete limitTree; 
	delete inputTree;
}

void Process::Run() {
	PInfo("LimitTreeBuilder::Process::Run",TString::Format("%s%s",name.Data(),syst.Data()));

	inputTree->SetBranchStatus("*",0);
	turnOnBranches(inputTree,selection);
	turnOnBranches(inputTree,weight);
	for (auto *x : vars) {
		turnOnBranches(inputTree,x->formula);
	}
	for (auto *x : formulae) {
		turnOnBranches(inputTree,x->formula);
	}

	// apply selection
	TTree *clonedTree = (TTree*)inputTree->CopyTree(selection.Data());

	// load inputs	
	for (auto *x : vars) 
		clonedTree->SetBranchAddress(x->formula,x->val);

	std::vector<TTreeFormula*> treeformulae;
	for (auto *x : formulae) {
		TTreeFormula *tf = new TTreeFormula(x->name.Data(),x->formula.Data(),clonedTree);
		tf->SetQuickLoad(true);
		tf->GetNdata();
		treeformulae.push_back(tf);
	}
	unsigned int nF = treeformulae.size();

	TTreeFormula fweight(TString::Format("w_%s",name.Data()).Data(),weight.Data(),clonedTree);
	fweight.SetQuickLoad(true);
	fweight.GetNdata();
	float weightval=0; 
	limitTree->Branch("weight",&weightval,"weight/F");

	// loop through and do stuff
	unsigned int nEntries = clonedTree->GetEntries(), iE=0;
	ProgressReporter pr("LimitTreeBuilder::Process::Run",&iE,&nEntries,10);
	for (iE=0; iE!=nEntries; ++iE) {
		pr.Report();
		clonedTree->GetEntry(iE);
		weightval = fweight.EvalInstance();
		for (unsigned int iF=0; iF!=nF; ++iF) {
			*(formulae[iF]->val) = treeformulae[iF]->EvalInstance();
		}
		limitTree->Fill();
	}

	for (auto tf : treeformulae)
		delete tf;
	delete clonedTree;
}

void LimitTreeBuilder::Output() {
	for (auto r : regions) {
		const char *rname = r->name.Data();
		for (auto p : r->GetProcesses()) {
			const char *pname = p->name.Data();
			const char *systname = p->syst.Data();
			fOut->WriteTObject(p->GetTree(),TString::Format("%s_%s%s",pname,rname,systname));
		}
	}
	fOut->Close();
	fOut=0;
}
