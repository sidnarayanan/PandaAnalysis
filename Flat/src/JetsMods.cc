#include "../interface/JetsMods.h"

using namespace pa;
using namespace std;
using namespace panda; 

void JetMod::do_readData(TString dirPath) {
  if (!analysis.rerunJES)
    return;

  TString jecV = "V4", jecReco = "23Sep2016"; 
  TString jecVFull = jecReco+jecV;
  std::vector<TString> eraGroups = {"BCD","EF","G","H"};

  ak4UncReader["MC"] = new JetCorrectionUncertainty(
       (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_Uncertainty_AK4PFPuppi.txt").Data()
    );
  for (auto e : eraGroups) {
    ak4UncReader["data"+e] = new JetCorrectionUncertainty(
         (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_Uncertainty_AK4PFPuppi.txt").Data()
      );
  }

  ak4JERReader = new JERReader(dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_SF_AK4PFPuppi.txt",
                               dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_PtResolution_AK4PFPuppi.txt");

  std::vector<JetCorrectorParameters> params = {
    JetCorrectorParameters(
      (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L1FastJet_AK4PFPuppi.txt").Data()),
    JetCorrectorParameters(
      (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L2Relative_AK4PFPuppi.txt").Data()),
    JetCorrectorParameters(
      (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L3Absolute_AK4PFPuppi.txt").Data()),
    JetCorrectorParameters(
      (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_L2L3Residual_AK4PFPuppi.txt").Data())
  };
  ak4ScaleReader["MC"] = new FactorizedJetCorrector(params);
  if (DEBUG>1) PDebug("JetMod::do_readData","Loaded JES for AK4 MC");
  for (auto e : eraGroups) {
    params = {
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L1FastJet_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L2Relative_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L3Absolute_AK4PFPuppi.txt").Data()),
      JetCorrectorParameters(
        (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_L2L3Residual_AK4PFPuppi.txt").Data())
    };
    ak4ScaleReader["data"+e] = new FactorizedJetCorrector(params);
    if (DEBUG>1) PDebug("JetMod::do_readData","Loaded JES for AK4 "+e);
  }

  if (DEBUG) PDebug("JetMod::do_readData","Loaded JES/R");
}

void JetMod::do_init(Registry& registry) 
{
  jesShifts = registry.access<std::vector<JESHandler>>("jesShifts");
}

void JetMod::setupJES()
{
  if (!analysis.rerunJES || (uncReaderAK4 != nullptr)) 
    return;
  if (analysis.isData) {
    TString thisEra = utils.eras.getEra(gt.runNumber);
    for (auto& iter : ak4UncReader) {
      if (!iter.first.Contains("data"))
        continue;
      if (iter.first.Contains(thisEra)) {
        uncReaderAK4 = ak4UncReader[iter.first];
        scaleReaderAK4 = ak4ScaleReader[iter.first];
        return;
      }
    }
  } else {
    uncReaderAK4 = ak4UncReader["MC"];
    scaleReaderAK4 = ak4ScaleReader["MC"];
  }
}

void JetMod::do_execute()
{
  setupJES();


}