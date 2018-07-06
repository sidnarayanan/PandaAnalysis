#include "PandaAnalysis/Utilities/interface/EnergyCorrelations.h"
#include "PandaAnalysis/Utilities/interface/CSVHelper.h"
#include "PandaAnalysis/Utilities/interface/RoccoR.h"
#include "PandaAnalysis/Utilities/interface/RelIso.h"
#include "PandaAnalysis/Utilities/interface/EtaPhiMap.h"
#include "PandaAnalysis/Utilities/interface/PackingHelperStandalone.h"
#include "PandaAnalysis/Utilities/interface/JetCorrector.h"
#include "PandaAnalysis/Utilities/interface/KinematicFit.h"

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace pa;
#pragma link C++ namespace kinfit;

#pragma link C++ class pa::ECFCalculator;
#pragma link C++ class pa::ECFCalculator::iterator;
#pragma link C++ class pa::RocRes;
#pragma link C++ class pa::RocOne;
#pragma link C++ class pa::RoccoR;
#pragma link C++ class pa::CSVHelper;
#pragma link C++ class pa::EtaPhiMap;
#pragma link C++ class pa::MiniIsoParams;
#pragma link C++ class pa::PackingHelperStandalone;
#pragma link C++ class pa::JetCorrector;
#pragma link C++ class kinfit::Fit;


#endif
