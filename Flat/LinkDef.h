#include "PandaAnalysis/Flat/interface/AnalyzerUtilities.h"
#include "PandaAnalysis/Flat/interface/GeneralTree.h"
#include "PandaAnalysis/Flat/interface/TagTree.h"
#include "PandaAnalysis/Flat/interface/JetCorrector.h"
#include "PandaAnalysis/Flat/interface/LimitTreeBuilder.h"
#include "PandaAnalysis/Flat/interface/PandaAnalyzer.h"
#include "PandaAnalysis/Flat/interface/TagAnalyzer.h"
#include "PandaAnalysis/Flat/interface/genericTree.h"


#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ enum PandaAnalyzer::PreselectionBit;
#pragma link C++ enum PandaAnalyzer::ProcessType;
#pragma link C++ enum PandaAnalyzer::TriggerBits;
#pragma link C++ enum TagAnalyzer::ProcessType;
#pragma link C++ enum GeneralTree::BTagShift;
#pragma link C++ enum GeneralTree::BTagJet;
#pragma link C++ enum GeneralTree::BTagTags;

#pragma link C++ class Analysis;
#pragma link C++ class LumiRange;
#pragma link C++ class TCorr;
#pragma link C++ class THCorr;
#pragma link C++ class TF1Corr;
#pragma link C++ class btagcand;
#pragma link C++ class JetCorrector;
#pragma link C++ class PandaAnalyzer;
#pragma link C++ class TagAnalyzer;
#pragma link C++ class GeneralTree;
#pragma link C++ class GeneralTree::ECFParams;
#pragma link C++ class GeneralTree::BTagParams;
#pragma link C++ class TagTree;
#pragma link C++ class TagTree::ECFParams;
#pragma link C++ class genericTree;
#pragma link C++ class LimitTreeBuilder;
#pragma link C++ class xformula;
#pragma link C++ class VariableMap;
#pragma link C++ class Process;
#pragma link C++ class Region;
#pragma link C++ class ParticleGridder;
#pragma link C++ class JetTree;

#endif
