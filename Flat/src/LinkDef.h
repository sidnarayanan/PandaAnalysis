#include "PandaAnalysis/Flat/interface/AnalyzerUtilities.h"
#include "PandaAnalysis/Flat/interface/GeneralTree.h"
#include "PandaAnalysis/Flat/interface/TagTree.h"
#include "PandaAnalysis/Flat/interface/LimitTreeBuilder.h"
#include "PandaAnalysis/Flat/interface/PandaAnalyzer.h"
#include "PandaAnalysis/Flat/interface/HRAnalyzer.h"
#include "PandaAnalysis/Flat/interface/TagAnalyzer.h"
#include "PandaAnalysis/Flat/interface/genericTree.h"
#include "PandaAnalysis/Flat/interface/Selection.h"
#include "PandaAnalysis/Flat/interface/L1Analyzer.h"


#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace pa;

#pragma link C++ enum pa::ProcessType;
#pragma link C++ enum pa::TriggerBits;
#pragma link C++ enum GeneralTree::BTagShift;
#pragma link C++ enum GeneralTree::BTagJet;
#pragma link C++ enum GeneralTree::BTagTags;

#pragma link C++ class pa::Selection;
#pragma link C++ class pa::LambdaSel;
#pragma link C++ class pa::TriggerSel;
#pragma link C++ class pa::GenBosonPtSel;
#pragma link C++ class pa::FatJetSel;
#pragma link C++ class pa::FatJet450Sel;
#pragma link C++ class pa::GenFatJetSel;
#pragma link C++ class pa::LeptonSel;
#pragma link C++ class pa::LeptonFakeSel;
#pragma link C++ class pa::RecoilSel;
#pragma link C++ class pa::MonotopSel;
#pragma link C++ class pa::MonohiggsSel;
#pragma link C++ class pa::VHbbSel;
#pragma link C++ class pa::Analysis;
#pragma link C++ class pa::LumiRange;
#pragma link C++ class pa::TCorr;
#pragma link C++ class pa::THCorr;
#pragma link C++ class pa::TF1Corr;
#pragma link C++ class pa::btagcand;
#pragma link C++ class pa::PandaAnalyzer;
#pragma link C++ class pa::HRAnalyzer;
#pragma link C++ class pa::TagAnalyzer;
#pragma link C++ class pa::GeneralTree;
#pragma link C++ class pa::GeneralTree::ECFParams;
#pragma link C++ class pa::GeneralTree::BTagParams;
#pragma link C++ class pa::TagTree;
#pragma link C++ class pa::TagTree::ECFParams;
#pragma link C++ class pa::genericTree;
#pragma link C++ class pa::LimitTreeBuilder;
#pragma link C++ class pa::xformula;
#pragma link C++ class pa::VariableMap;
#pragma link C++ class pa::Process;
#pragma link C++ class pa::Region;
#pragma link C++ class pa::ParticleGridder;
#pragma link C++ class pa::JetTree;
#pragma link C++ class pa::L1Analyzer; 

#endif
