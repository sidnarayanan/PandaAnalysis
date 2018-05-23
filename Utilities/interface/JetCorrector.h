#ifndef PANDACORE_TOOLS_JetCorrector
#define PANDACORE_TOOLS_JetCorrector

#include "PandaCore/Tools/interface/Common.h"
#include "PandaCore/Tools/interface/DataTools.h"
#include <map>
#include <string>
#include "TString.h"
#include "PandaTree/Objects/interface/Jet.h"
#include "PandaTree/Objects/interface/Met.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

/**
 * \brief Corrects a jet collection and optionally propagates to Met 
 *
 */
namespace pa {
	class JetCorrector
	{
	public:
		JetCorrector();
		~JetCorrector();

		void RunCorrection(bool isData, float rho, panda::JetCollection *injets_, panda::Met *rawmet_=0, int runNumber = 0);
		panda::JetCollection *GetCorrectedJets();
		panda::Met *GetCorrectedMet();

		void SetMCCorrector(TString fpath);
		void SetDataCorrector(TString fpath, TString iov = "all");

	private:
			FactorizedJetCorrector *mMCJetCorrector;	 
			std::map<TString,FactorizedJetCorrector *> mDataJetCorrectors;	// map from era to corrector

			panda::JetCollection *outjets = 0;
			panda::Met *outmet = 0;

			EraHandler *era = 0;	
	};
}
#endif
