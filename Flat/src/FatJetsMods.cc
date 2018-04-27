#include "../interface/FatJetsMods.h"

using namespace pa;
using namespace std;
using namespace panda; 

void FatJetMod::do_readData(TString dirPath) {
  if (!analysis.rerunJES)
    return;

  TString jecV = "V4", jecReco = "23Sep2016"; 
  TString jecVFull = jecReco+jecV;
  std::vector<TString> eraGroups = {"BCD","EF","G","H"};

  ak8UncReader["MC"] = new JetCorrectionUncertainty(
       (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecVFull+"_MC_Uncertainty_AK8PFPuppi.txt").Data()
    );
  for (auto e : eraGroups) {
    utils.ak8UncReader["data"+e] = new JetCorrectionUncertainty(
         (dirPath+"/jec/"+jecVFull+"/Summer16_"+jecReco+e+jecV+"_DATA_Uncertainty_AK8PFPuppi.txt").Data()
      );
  }

  ak8JERReader = new JERReader(dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_SF_AK8PFPuppi.txt",
                               dirPath+"/jec/25nsV10/Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt");

}


void FatJetMod::setupJES()
{
  if (!analysis.rerunJES || (uncReader != nullptr)) 
    return;
  if (analysis.isData) {
    TString thisEra = utils.eras.getEra(gt.runNumber);
    for (auto& iter : ak8UncReader) {
      if (!iter.first.Contains("data"))
        continue;
      if (iter.first.Contains(thisEra)) {
        uncReader = ak8UncReader[iter.first];
        return;
      }
    }
  } else {
    uncReader = ak8UncReader["MC"];
  }
}

void FatJetMod::do_execute()
{
  setupJES();

  gt.nFatjet=0;
  int fatjet_counter=-1;
  for (auto& fj : *fatjets) {
    ++fatjet_counter;
    float pt = fj.pt();
    float rawpt = fj.rawPt;
    float eta = fj.eta();
    float mass = fj.m();
    float ptcut = 200;
    if (analysis.deep)
      ptcut = 400;
    
    float bestPt = pt;
    if (analysis.rerunJES) {
      bestPt = TMath::Max(bestPt, fj.ptCorrUp);
      bestPt = TMath::Max(bestPt, fj.ptCorrDown);
    }
    if (bestPt<ptcut || fabs(eta)>2.4 || !fj.monojet)
      continue;

    float phi = fj.phi();
    if (isMatched(matchLeps,FATJETMATCHDR2,eta,phi) || 
        isMatched(matchPhos,FATJETMATCHDR2,eta,phi)) {
      continue;
    }

    gt.nFatjet++;
    if (gt.nFatjet==1) {
      fj1 = &fj;
      gt.fj1IsClean = fatjet_counter==0 ? 1 : 0;
      gt.fjPt = pt;
      gt.fjEta = eta;
      gt.fjPhi = phi;
      gt.fjM = mass;
      gt.fjMSD = fj.mSD;
      gt.fjRawPt = rawpt;

      // mSD correction
      float corrweight=1.;
      corrweight = getMSDCorr(pt,eta);
      gt.fjMSD_corr = corrweight*gt.fjMSD;

      // now we do substructure
      gt.fjTau32 = clean(fj.tau3/fj.tau2);
      gt.fjTau32SD = clean(fj.tau3SD/fj.tau2SD);
      gt.fjTau21 = clean(fj.tau2/fj.tau1);
      gt.fjTau21SD = clean(fj.tau2SD/fj.tau1SD);

      for (auto ibeta : ibetas) {
        for (auto N : Ns) {
          for (auto order : orders) {
            GeneralTree::ECFParams p;
            p.order = order; p.N = N; p.ibeta = ibeta;
            if (gt.fjIsClean || true)
              gt.fjECFNs[p] = fj.get_ecf(order,N,ibeta);
            else
              gt.fjECFNs[p] = fj.get_ecf(order,N,ibeta);
          }
        }
      } //loop over betas
      gt.fjHTTMass = fj.htt_mass;
      gt.fjHTTFRec = fj.htt_frec;

      std::vector<panda::MicroJet const*> subjets;
      for (int iS(0); iS != fj.subjets.size(); ++iS)
        subjets.push_back(&fj.subjets.objAt(iS));

      auto csvsort = [](panda::MicroJet const* j1, panda::MicroJet const* j2) -> bool {
              return j1->csv > j2->csv;
            };

      std::sort(subjets.begin(),subjets.end(),csvsort);
      if (subjets.size()>0) {
        gt.fjMaxCSV = subjets.at(0)->csv;
        gt.fjMinCSV = subjets.back()->csv;
        if (subjets.size()>1) {
          gt.fjSubMaxCSV = subjets.at(1)->csv;
        }
      }

      gt.fjDoubleCSV = fj.double_sub;
      if (analysis.monoh) {
        for (int iSJ=0; iSJ!=fj.subjets.size(); ++iSJ) {
          auto& subjet = fj.subjets.objAt(iSJ);
          gt.fjsjPt[iSJ]=subjet.pt();
          gt.fjsjEta[iSJ]=subjet.eta();
          gt.fjsjPhi[iSJ]=subjet.phi();
          gt.fjsjM[iSJ]=subjet.m();
          gt.fjsjCSV[iSJ]=subjet.csv;
          gt.fjsjQGL[iSJ]=subjet.qgl;
        }
      }
    }
    if (!isData && fj.matchedGenJet.isValid())
      gt.fjGenNumB = fj.matchedGenJet.get()->numB;
    else 
      gt.fjGenNumB = 0;
  }

  recluster->execute();
}


float FatJetMod::getMSDCorr(float puppipt, float puppieta) 
{

  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;

  genCorr = utils.puppisd_corrGEN->Eval( puppipt );
  if ( fabs(puppieta) <= 1.3 ){
    recoCorr = utils.puppisd_corrRECO_cen->Eval( puppipt );
  }
  else {
    recoCorr = utils.puppisd_corrRECO_for->Eval( puppipt );
  }
  totalWeight = genCorr * recoCorr;

  return totalWeight;
}

void FatJetReclusterMod::do_execute()
{
  if ((*fj1) == nullptr)
    return;

  VPseudoJet particles = convertPFCands(event.pfCandidates,analysis.puppiJets,0);
  ClusterSequenceArea seq(particles,*jetDef,*(utils.areaDef));
  VPseudoJet allJets(seq.inclusive_jets(0.));
  fastjet::PseudoJet *pj1=0;
  double minDR2 = 999;
  for (auto &jet : allJets) {
    double dr2 = DeltaR2(jet.eta(),jet.phi_std(),(*fj1)->eta(),(*fj1)->phi());
    if (dr2<minDR2) {
      minDR2 = dr2;
      pj1 = &jet;
    }
  }
  if (pj1) {
    VPseudoJet constituents = fastjet::sorted_by_pt(pj1->constituents());

    double eTot=0, eTrunc=0;
    for (unsigned iC=0; iC!=constituents.size(); ++iC) {
      double e = constituents.at(iC).E();
      eTot += e;
      if (iC<100)
        eTrunc += e;
    }


    fastjet::PseudoJet sdJet = (*softDrop)(*pj1);
    VPseudoJet sdConstituents = fastjet::sorted_by_pt(sdJet.constituents());
    eTot=0; eTrunc=0;
    for (unsigned iC=0; iC!=sdConstituents.size(); ++iC) {
      double e = sdConstituents.at(iC).E();
      eTot += e;
      if (iC<100)
        eTrunc += e;
    }
}
