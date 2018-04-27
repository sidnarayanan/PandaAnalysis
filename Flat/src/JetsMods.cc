#include "../interface/JetsMods.h"

using namespace pa;
using namespace std;
using namespace panda; 

float centralOnly(float x, float aeta, float def = -1) 
{
  return  aeta < 2.4 ? x : -1;
}

bool csvLoose(float csv) 
{
  return csv > 0.5426;
}

bool csvMed(float csv) 
{
  return csv > 0.8484;
}

JetWrapper shiftJet(const Jet& jet, shiftjes shift, bool smear=false) 
{
  float pt;
  switch (shift) {
    case shiftjes::kNominal:
      pt = jet.pt();
      break;
    case shiftjes::kJESUp:
      pt = jet.ptCorrUp;
      break;
    case shiftjes::kJESDown:
      pt = jet.ptCorrDown;
      break;
    default:
      PError("shiftJet", "Unknown JES type!");
      exit(1);
  }
  if (smear) {
    pt *= jet.ptSmear / jet.pt();
  }
  return JetWrapper(pt, jet);
}


void JetMod::do_readData(TString dirPath) 
{
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

void JetMod::varyJES()
{
  JESLOOP {
    auto& jets = (*jesShifts)[shift];
    jets.reserve(ak4jets->size());
    for (auto &j : *ak4jets) {
      jets.all.push_back(shiftJet(j, i2jes(shift), analysis.hbb));
    }
    if (shift != jes2i(shiftjes::kNominal)) {
      std::sort(jets.all.begin(), jets.all.end(),
                [](JetWrapper x, JetWrapper y) { return x.pt > y.pt; });
    }
  }
}

void JetMod::do_execute()
{
  setupJES();

  varyJES();

  float maxJetEta = analysis.vbf ? 4.7 : 4.5;
  int nJetDPhi = analysis.vbf ? 4 : 5;
  float minMinJetPt = min(cfg.minJetPt, cfg.minBJetPt);

  TLorentzVector vBarrelJets;

  JESLOOP {
    JESHandler& jets = (*jesShifts)[shift];
    for (auto& jw : jets.all) {
      auto& jet = jw.get_base();
      float aeta = abs(jet.eta());
      float pt = jw.pt;
      if (aeta > maxJetEta || pt < minMinJetPt) 
        continue;
      if (IsMatched(matchLeps,0.16,jet.eta(),jet.phi()))
        continue;
      if (!analysis.hbb && IsMatched(matchPhos,0.16,jet.eta(),jet.phi()))
        continue;
      if ((analysis.vbf || analysis.hbb) && !jet.loose)
        continue;

      if (analysis.jetFlavorPartons)
        jetPartonFlavor(jet, jw.flavor, jw.genpt);
      else if (analysis.jetFlavorJets)
        jetClusteredFlavor(jet, jw.flavor, jw.genpt);

      float csv = centralOnly(jet.csv, aeta);
      float cmva = centralOnly(jet.cmva, aeta);

      if (pt > minBJetPt && aeta < 2.4) { // b jets
        if (csvLoose(csv)) { 
          ++(gt.jetNBtags[shift]);
          if (csvMed(csv)) {
            ++(gt.jetNMBtags[shift]);
          }
        }

        jets.bcand.push_back(&jw); 
      } 


      if (pt > minJetPt) {
        if ((analysis.hbb || analysis.monoh) && jets.cleaned.size() >= NJET)
          continue;

        jets.cleaned.push_back(&jw);
        
        if (jets.cleaned.size() < 3) {      
          if (GetCorr(cBadECALJets, jet.eta(), jet.phi()) > 0)
            gt.badECALFilter = 0;
        }

        if (analysis.fatjet) {
          IsoJet(jw, jets);
          if (jw.iso && csvLoose(csv))
            ++(gt.isojetNBtags[shift]);
        }

        if (aeta < 2.4) {
          jets.central.push_back(&jw);
          
          int njet = jets.central.size();
          gt.nJet[shift] = njet;
          if (shift == jes2i(shiftjes::kNominal)) {
            if (njet < 2) {
              gt.jetPt[njet] = pt;
              gt.jetEta[njet] = jet.eta();
              gt.jetPhi[njet] = jet.phi();
              gt.jetCSV[njet] = csv;
              gt.jetIsTight[njet] = jet.monojet ? 1 : 0;
              gt.jetFlav[njet] = jw.flavor;
              gt.jetIsIso[njet] = jw.iso ? 1 : 0;
            }
          }
        }

        if (shift == jes2i(shiftjes::kNominal) && aeta < 3.0) {
          gt.barrelHT += jet.pt();
          vBarrelJets += jet.p4();
        }

        int njet = jets.cleaned.size();
        if (njet < 2 || ((analysis.hbb || analysis.monoh) && njet < NJET)) {
          gt.jotPt[njet][shift] = pt;
          gt.jotEta[njet][shift] = jet.eta();
          gt.jotPhi[njet][shift] = jet.phi();
          gt.jotE[njet][shift] = jet.e() * jw.scale();
          gt.jotCSV[njet][shift] = csv;
          gt.jotCMVA[njet][shift] = cmva;
          gt.jotVBFID[njet][shift] = (aeta < 2.4) ? (jet.monojet ? 1 : 0) : 1;

          if (analysis.bjetRegression) {
            JetBRegressionInfo(jet, njet, shift);
          }
        }

        TLorentzVector vJet; vJet.SetPtEtaPhiM(pt, jet.eta(), jet.phi(), jet.m());
        if (njet < nJetDPhi) {
          gt.dphipuppimet[shift] = min(fabs(vJet.DeltaPhi(jets.vpuppiMET)), (double)gt.dphipuppimet[shift]);
          gt.dphipfmet[shift]    = min(fabs(vJet.DeltaPhi(jets.vpfMET)),    (double)gt.dphipfmet[shift]);
          if (analysis.recoil) {
            gt.dphipuppiUA[shift] = min(fabs(vJet.DeltaPhi(jets.vpuppiUA)), (double)gt.dphipuppiUA[shift]);
            gt.dphipuppiUW[shift] = min(fabs(vJet.DeltaPhi(jets.vpuppiUW)), (double)gt.dphipuppiUW[shift]);
            gt.dphipuppiUZ[shift] = min(fabs(vJet.DeltaPhi(jets.vpuppiUZ)), (double)gt.dphipuppiUZ[shift]);
            gt.dphipfUA[shift] = min(fabs(vJet.DeltaPhi(jets.vpfUA)), (double)gt.dphipfUA[shift]);
            gt.dphipfUW[shift] = min(fabs(vJet.DeltaPhi(jets.vpfUW)), (double)gt.dphipfUW[shift]);
            gt.dphipfUZ[shift] = min(fabs(vJet.DeltaPhi(jets.vpfUZ)), (double)gt.dphipfUZ[shift]);
          }
        }
      }
    }

    gt.nJot[shift] = jets.cleaned.size();
    switch (gt.whichRecoil) {
      case 0: // MET
        gt.dphipuppiU[shift] = gt.dphipuppimet[shift];
        gt.dphipfU[shift] = gt.dphipfmet[shift];
        break;
      case -1: // photon
        gt.dphipuppiU[shift] = gt.dphipuppiUA[shift];
        gt.dphipfU[shift] = gt.dphipfUA[shift];
        break;
      case 1:
        gt.dphipuppiU[shift] = gt.dphipuppiUW[shift];
        gt.dphipfU[shift] = gt.dphipfUW[shift];
        break;
      case 2:
        gt.dphipuppiU[shift] = gt.dphipuppiUZ[shift];
        gt.dphipfU[shift] = gt.dphipfUZ[shift];
        break;
      default: // c'est impossible !
        break;
    }

    if (analysis.vbf) 
      JetVBFSystem(shift);

    if (analysis.monoh || analysis.hbb)
      JetHbbReco(shift);

  } // shift loop 

  gt.barrelHTMiss = vBarrelJets.Pt();
  
}
