#include "../interface/JetsOps.h"
#include "PandaAnalysis/Utilities/interface/Helicity.h"
#include "PandaAnalysis/Utilities/interface/NeutrinoSolver.h"
#include "TSystem.h"

using namespace pa;
using namespace std;
using namespace panda;
namespace fj = fastjet;
using JECParams = JetCorrectorParameters;

inline float centralOnly(float x, float aeta, float def = -1)
{
  return  aeta < 2.4 ? x : -1;
}

JetWrapper BaseJetOp::shiftJet(const Jet& jet, shiftjes shift, bool smear)
{
  float pt = jet.pt();
  if (smear) {
    if (recalcJER) {
      double smearFac=1, smearFacUp=1, smearFacDown=1;
      jer->getStochasticSmear(pt,jet.eta(),event.rho,smearFac,smearFacUp,smearFacDown);
      pt *= smearFac;
    } else {
      pt = jet.ptSmear;
    }
  }
  if (shift != shiftjes::kNominal) {
    int ishift = jes2i(shift);
    bool isUp = !(ishift % 2 == 0);
    if (analysis.rerunJES) {
      (*scaleUnc)[ishift]->setJetPt(pt);
      (*scaleUnc)[ishift]->setJetEta(jet.eta());
      double relShift = (*scaleUnc)[ishift]->getUncertainty(isUp);
      if (!isUp)
        relShift = -relShift; 
      pt *= (1 + relShift);
    } else {
      pt = (isUp ? jet.ptCorrUp :  jet.ptCorrDown) * pt / jet.pt();
    }
  }
  return JetWrapper(pt, jet);
}


void BaseJetOp::do_readData(TString dirPath)
{
  if (recalcJER) {
    jer.reset(new JERReader(dirPath+"/jec/"+jerV+"/"+jerV+"_MC_SF_"+jetType+".txt",
                            dirPath+"/jec/"+jerV+"/"+jerV+"_MC_PtResolution_"+jetType+".txt"));
  }

  if (!analysis.rerunJES)
    return;

  TString jecVFull = jecReco+spacer+jecV;

  TString basePath = dirPath+"/jec/"+jecVFull+"/"+campaign+"_"+jecVFull;
  vector<JECParams> params = {
    JECParams((basePath+"_MC_L1FastJet_"+jetType+".txt").Data()),
    JECParams((basePath+"_MC_L2Relative_"+jetType+".txt").Data()),
    JECParams((basePath+"_MC_L3Absolute_"+jetType+".txt").Data()),
    JECParams((basePath+"_MC_L2L3Residual_"+jetType+".txt").Data())
  };
  scales["MC"].reset(new FactorizedJetCorrector(params));
  for (auto e : eraGroups) {
    basePath = dirPath+"/jec/"+jecVFull+"/"+campaign+"_"+jecReco+e+spacer+jecV;
    params = {
      JECParams((basePath+"_DATA_L1FastJet_"+jetType+".txt").Data()),
      JECParams((basePath+"_DATA_L2Relative_"+jetType+".txt").Data()),
      JECParams((basePath+"_DATA_L3Absolute_"+jetType+".txt").Data()),
      JECParams((basePath+"_DATA_L2L3Residual_"+jetType+".txt").Data())
    };
    scales["data"+e].reset(new FactorizedJetCorrector(params));
  }


  basePath = dirPath+"/jec/"+jecVFull+"/"+campaign+"_";
  setScaleUnc("MC", (basePath+jecVFull+"_MC_UncertaintySources_"+jetType+".txt").Data());
  for (auto e : eraGroups) {
    setScaleUnc("data"+e, (basePath+jecReco+e+spacer+jecV+"_DATA_UncertaintySources_"+jetType+".txt").Data());
  }
}


void BaseJetOp::setScaleUnc(TString tag, TString path)
{
  scaleUncs[tag] = std::vector<std::shared_ptr<JetCorrectionUncertainty>>(jes2i(shiftjes::N),nullptr);
  JESLOOP {
    if (shift % 2 == 0)
      continue;
    TString shiftName = jesName(i2jes(shift));
    shiftName.ReplaceAll("JES","");
    shiftName = shiftName(0,shiftName.Length()-2);
    scaleUncs[tag][shift] = make_shared<JetCorrectionUncertainty>(JECParams(path.Data(), shiftName.Data()));
    scaleUncs[tag][shift+1] = scaleUncs[tag][shift]; // Down = Up
  }
}

void JetOp::setupJES()
{
  if (!analysis.rerunJES || (scaleUnc != nullptr))
    return;
  if (analysis.isData) {
    TString thisEra = utils.eras->getEra(gt.runNumber);
    for (auto& iter : scaleUncs) {
      if (!iter.first.Contains("data"))
        continue;
      if (iter.first.Contains(thisEra)) {
        scaleUnc = &(scaleUncs[iter.first]);
        scale = scales[iter.first].get();
        return;
      }
    }
  } else {
    scaleUnc = &(scaleUncs["MC"]);
    scale = scales["MC"].get();
  }
}


void JetOp::varyJES()
{
  JESLOOP {
    auto& jets = (*jesShifts)[shift];
    jets.reserve(ak4Jets->size());
    for (auto &j : *ak4Jets) {
      jets.all.push_back(shiftJet(j, i2jes(shift), analysis.hbb && !analysis.isData));
    }
  }
  for (size_t iJ = 0; iJ != (*jesShifts)[0].all.size(); ++iJ) {
    auto* nominal = &((*jesShifts)[0].all[iJ]);
    nominal->maxpt = 0;
    JESLOOP {
      auto* jet = &((*jesShifts)[shift].all[iJ]);
      jet->nominal = nominal;
      if (jet->pt > nominal->maxpt)
        nominal->maxpt = jet->pt;
    }
  }
}


void JetOp::do_execute()
{
  setupJES();

  varyJES();

  float maxJetEta = analysis.vbf ? 4.7 : 4.5;
  int nJetDPhi = analysis.vbf ? 4 : 5;
  float minMinJetPt = min(cfg.minJetPt, cfg.minBJetPt);

  TLorentzVector vBarrelJets;

  JESLOOP {
    bool isNominal = (shift == jes2i(shiftjes::kNominal));
    bool metShift = (i2jes(shift) <= shiftjes::kJESTotalDown);
    JESHandler& jets = (*jesShifts)[shift];
    (*currentJES) = &jets;
    for (auto& jw : jets.all) {
      (*currentJet) = &jw;
      auto& jet = jw.get_base();
      float aeta = abs(jet.eta());
      float pt = jw.pt;

      if (analysis.year == 2016 || analysis.year == 2017) {
        if (isNominal && !isMatched(matchVeryLoosePhos.get(),0.16,jet.eta(),jet.phi())) {
          // prefiring weights
          gt.sf_l1Prefire *= (1.0 - utils.getCorr(cL1PreFiring,jet.eta(),pt));
        }
      }

      if (aeta > maxJetEta || jw.nominal->maxpt < minMinJetPt)
        continue;
      
      if (isNominal) { // perform cleaning only on the nominal jet and save flags
        if (isMatched(matchLeps.get(),0.16,jet.eta(),jet.phi()))
          jw.isLep = true;
        if (!analysis.hbb && isMatched(matchPhos.get(),0.16,jet.eta(),jet.phi()))
          jw.isPho = true;
        if (analysis.hbb && jet.puid < utils.getCorr(cJetLoosePUID, aeta,min(39.99f,pt)))
          jw.isPileupJet = true;
      }
      if (jw.nominal->isLep || jw.nominal->isPho || jw.nominal->isPileupJet)
        continue;
      if ((analysis.vbf || analysis.hbb) && !jet.loose)
        continue;

      float csv = centralOnly( (analysis.useDeepCSV? jet.deepCSVb + jet.deepCSVbb : jet.csv), aeta);
      float cmva = centralOnly(jet.cmva, aeta);

      if (pt > cfg.minBJetPt && aeta < 2.4) { // b jets
        if (csvLoose(csv)) {
          ++(gt.jetNBtags[shift]);
          if (csvMed(csv)) {
            ++(gt.jetNMBtags[shift]);
          }
        }

        if (isNominal)
          jets.bcand.push_back(&jw);
      }


      if (jw.nominal->maxpt > cfg.minJetPt) {
        // for H->bb, don't consider any jet past NJETSAVED, 
        // for other analyses, consider them, just don't save them
        if ((analysis.hbb || analysis.monoh) && (int)jets.cleaned.size() >= cfg.NJETSAVED)
          continue;

        jets.cleaned.push_back(&jw);

        if (jets.cleaned.size() < 3) {
          if (utils.getCorr(cBadECALJets, jet.eta(), jet.phi()) > 0)
            gt.badECALFilter = 0;
        }

        if (isNominal && analysis.fatjet) {
          isojet->execute();
          if (jw.iso && csvLoose(csv)) {
            ++(gt.isojetNBtags[shift]);
            if (csvMed(csv))
              ++(gt.isojetNMBtags[shift]);
          }
        }

        if (aeta < 2.4) {
          jets.central.push_back(&jw);

          int njet = jets.central.size() - 1;
          if (pt > cfg.minJetPt)
            gt.nJet[shift]++;
          if (isNominal) {
            if (njet < 2) {
              jw.central_idx = njet; 
              gt.jetPt[njet] = pt;
              gt.jetEta[njet] = jet.eta();
              gt.jetPhi[njet] = jet.phi();
              gt.jetCSV[njet] = csv;
              gt.jetIsTight[njet] = jet.monojet ? 1 : 0;
              gt.jetIsIso[njet] = jw.iso ? 1 : 0;
            }
          }
        }

        if (isNominal && aeta < 3.0) {
          gt.barrelHT += jet.pt();
          vBarrelJets += jet.p4();
        }

        int njet = jets.cleaned.size() - 1;
        if (njet < 2 || ((analysis.hbb || analysis.monoh) && njet < cfg.NJETSAVED)) {
          jw.cleaned_idx = njet; 
          gt.jotPt[shift][njet] = pt;
          if (isNominal) {
            if (!analysis.hbb && jet.matchedGenJet.isValid())
              gt.jotGenPt[njet] = jet.matchedGenJet->pt(); 
            gt.jotSmear[njet] = jw.pt / jet.pt(); // smeared / nominal
            gt.jotEta[njet] = jet.eta();
            gt.jotPhi[njet] = jet.phi();
            gt.jotM[njet] = jet.m();
            gt.jotCSV[njet] = csv;
            gt.jotFlav[njet] = jw.flavor;
            gt.jotCMVA[njet] = cmva;
            gt.jotVBFID[njet] = (aeta < 2.4) ? (jet.monojet ? 1 : 0) : 1;
            gt.jotIso[njet] = jw.iso ? 1 : 0; 

            bjetreg->execute();
          }
        }
        if (pt > cfg.minJetPt)
          gt.nJot[shift]++;

        TLorentzVector vJet; vJet.SetPtEtaPhiM(pt, jet.eta(), jet.phi(), jet.m());
        if (metShift && njet < nJetDPhi && pt > cfg.minJetPt) { // only do this for fully-correlated shifts
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

    gt.nJotMax = max(gt.nJot[shift], gt.nJotMax);
    if (metShift) {
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
    }

    // dijet system
    if (metShift) {
      jets.sort();
      vbf->execute();
    }
    hbb->execute();

  } // shift loop
  gt.barrelHTMiss = vBarrelJets.Pt();

}


void JetFlavorOp::partonFlavor(JetWrapper& jw)
{
  auto& jet = jw.get_base();
  jw.flavor=0; jw.genpt=0;
  for (auto* genptr : *genP) {
    auto& gen = pToGRef(genptr);
    int apdgid = abs(gen.pdgid);
    if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
      continue;
    double dr2 = DeltaR2(jet.eta(),jet.phi(),gen.eta(),gen.phi());
    if (dr2<0.09) {
      jw.genpt = gen.pt();
      if (apdgid==4 || apdgid==5) {
        jw.flavor=apdgid;
        return;
      } else {
        jw.flavor=0;
      }
    }
  }
}


void JetFlavorOp::clusteredFlavor(JetWrapper& jw)
{
  auto& jet = jw.get_base();
  for (auto &gen : event.ak4GenJets) {
    if (DeltaR2(gen.eta(), gen.phi(), jet.eta(), jet.phi()) < 0.09) {
      int apdgid = abs(gen.pdgid);
      jw.genpt = gen.pt();
      if (apdgid == 4 || apdgid == 5) {
        jw.flavor = apdgid;
        return;
      } else {
        jw.flavor = 0;
      }
    }
  }
}


void JetFlavorOp::do_execute()
{
  for (auto& jw : (*jesShifts)[0].all) {
    if (analysis.jetFlavorPartons)
      partonFlavor(jw);
    else
      clusteredFlavor(jw);
    if (jw.cleaned_idx >= 0) {
      gt.jotFlav[jw.cleaned_idx] = jw.flavor; 
      if (jw.central_idx >= 0) 
        gt.jetFlav[jw.central_idx] = jw.flavor;
    }
  }
}


void IsoJetOp::do_execute()
{
  auto& jw = **currentJet;
  auto& jets = **currentJES;
  auto& jet = jw.get_base();
  float maxIsoEta = analysis.monoh ? 4.5 : 2.5;
  bool isIsoJet = (
        gt.nFatJet == 0 ||
        (fabs(jet.eta()) < maxIsoEta &&
         DeltaR2(gt.fjEta[0],gt.fjPhi[0],jet.eta(),jet.phi()) > cfg.FATJETMATCHDR2)
      );

  jw.iso = isIsoJet;

  if (isIsoJet)
    jets.iso.push_back(&jw);
}


const vector<double> BJetRegOp::Energies::dr2_bins {
  pow(0.05, 2), pow(0.1, 2), pow(0.2, 2), pow(0.3, 2), pow(0.4, 2)
};

void BJetRegOp::do_execute()
{
  auto& jw = **currentJet;
  auto& jet = jw.get_base();
  auto& jets = **currentJES;

  int N = jets.cleaned.size() - 1;

  energies.clear();

  TLorentzVector vjet(jet.p4());
  TLorentzVector vraw(vjet); vraw *= jet.rawPt / jet.pt(); 

  gt.jotEMF[N] = jet.cef + jet.nef;
  gt.jotHF[N] = jet.chf + jet.nhf;
  gt.jotCEF[N] = jet.cef;
  gt.jotNEF[N] = jet.nef;
  gt.jotCHF[N] = jet.chf;
  gt.jotNHF[N] = jet.nhf;
  gt.jotRawPt[N] = vraw.Pt();
  gt.jotRawMt[N] = vraw.Mt();
  gt.jotRawEt[N] = vraw.Et();
  gt.jotRawM[N] = vraw.M();
  gt.jotRawE[N] = vraw.E();
  gt.jotRho[N] = event.rho;
  gt.jotArea[N] = jet.area; 
  energies.jet_e = gt.jotRawE[N];

  float sumpt{0}, sumpt2{0};
  int leadingLepPdgId = 0;
  for (const auto& pf : jet.constituents) {
    if (!pf.isValid()) {
      // not sure why this is happening, but catch it 
      logger.warning("BJetRegOp::do_execute",Form("Cannot access PF at idx %i out of %i", pf.idx(), event.pfCandidates.size()));
      continue; 
    }
    TLorentzVector v(pf->p4());
    float dr2 = DeltaR2(v.Eta(), v.Phi(), vraw.Eta(), vraw.Phi());
    float pt = pf->pt();
    sumpt += pt; sumpt2 += pow(pt, 2);
    if (pt > 0.3)
      gt.jotNPt03[N]++;
    int pdgid = abs(pf->pdgId());
    if (pf->q() != 0) {
      gt.jotChTrk1Pt[N] = max(pt, gt.jotChTrk1Pt[N]);
      if (pdgid == 11 || pdgid == 13) {
        gt.jotNLep[N]++;
        if (pt > gt.jotLep1Pt[N]) {
          gt.jotLep1Pt[N] = pt;
          gt.jotLep1Eta[N] = pf->eta();
          gt.jotLep1Phi[N] = pf->phi();
          gt.jotLep1PtRel[N] = v.Perp(vjet.Vect()); 
          gt.jotLep1PtRelRaw[N] = v.Perp(vraw.Vect()); 
          gt.jotLep1PtRelRawInv[N] = vraw.Perp(v.Vect());
          gt.jotLep1DeltaR[N] = sqrt(dr2);
          leadingLepPdgId = pdgid;
          gt.jotLep1IsOther[N] = 0; 
        }
      }
    }
    gt.jotTrk1Pt[N] = max(pt, gt.jotTrk1Pt[N]);

    if (leadingLepPdgId == 11) 
      gt.jotLep1IsEle[N] = 1; 
    else if (leadingLepPdgId == 13)
      gt.jotLep1IsMu[N] = 1;  // if neither, then the default is IsOther == 1

    unsigned bin = lower_bound(Energies::dr2_bins.begin(), Energies::dr2_bins.end(), dr2)
                   - Energies::dr2_bins.begin();
    if (bin < Energies::dr2_bins.size()) {
      Energies::pftype pf_type{Energies::pne};
      if (pdgid == 22 || pdgid == 11)
        pf_type = Energies::pem;
      else if (pdgid == 13)
        pf_type = Energies::pmu;
      else if (pf->q() != 0)
        pf_type = Energies::pch;
      energies.pf[pf_type][bin].push_back(v);
    }
  }

  gt.jotPtD[N] = sumpt > 0 ? sqrt(sumpt2)/sumpt : 0;
  for (unsigned b = 0; b != Energies::dr2_bins.size(); ++b) {
    gt.jotEMRing[b][N] = energies.get_e(b, Energies::pem);
    gt.jotChRing[b][N] = energies.get_e(b, Energies::pch);
    gt.jotMuRing[b][N] = energies.get_e(b, Energies::pmu);
    gt.jotNeRing[b][N] = energies.get_e(b, Energies::pne);
  }
  /*
  for (int pf_type = Energies::pem; pf_type != (int)Energies::pN; ++pf_type) {
    static std::array<float, static_cast<long unsigned>(shiftjetrings::N)> moments;
    // eta first
    energies.get_moments(pf_type, &TLorentzVector::Eta, moments, jet.eta());
    auto eta_array = [&] () {
      switch (pf_type) {
        case Energies::pem: return &GeneralTree::jotEMEta;
        case Energies::pch: return &GeneralTree::jotChEta;
        case Energies::pmu: return &GeneralTree::jotMuEta;
        default:            return &GeneralTree::jotNeEta;
      }
    } ();
    for (int i = 0; i != static_cast<int>(shiftjetrings::N); ++i)
      (gt.*eta_array)[i][N] = moments[i];

    // phi second
    energies.get_moments(pf_type, &TLorentzVector::Phi, moments, jet.phi());
    auto phi_array = [&] () {
      switch (pf_type) {
        case Energies::pem: return &GeneralTree::jotEMPhi;
        case Energies::pch: return &GeneralTree::jotChPhi;
        case Energies::pmu: return &GeneralTree::jotMuPhi;
        default:            return &GeneralTree::jotNePhi;
      }
    } ();
    for (int i = 0; i != static_cast<int>(shiftjetrings::N); ++i)
      (gt.*phi_array)[i][N] = moments[i];

    auto dr_array = [&] () {
      switch (pf_type) {
        case Energies::pem: return &GeneralTree::jotEMDR;
        case Energies::pch: return &GeneralTree::jotChDR;
        case Energies::pmu: return &GeneralTree::jotMuDR;
        default:            return &GeneralTree::jotNeDR;
      }
    } ();
    for (int i = 0; i != static_cast<int>(shiftjetrings::N); ++i)
      (gt.*dr_array)[i][N] = sqrt(pow((gt.*eta_array)[i][N],2) + pow((gt.*phi_array)[i][N],2));
  }
  */

  auto& vert = jet.secondaryVertex;
  if (vert.isValid()) {
    gt.jotVtxPt[N] = vert->pt();
    gt.jotVtxMass[N] = vert->m();
    gt.jotVtx3DVal[N] = vert->vtx3DVal;
    gt.jotVtx3DErr[N] = vert->vtx3DeVal;
    gt.jotVtxNtrk[N] = vert->ntrk;
  }
}

void VBFSystemOp::do_execute()
{
  auto& jets = **currentJES;

  int shift = jets.shift_idx;

  unsigned idx0=0, idx1=1;
  if (analysis.hbb) {
    if (analysis.fatjet && fjPtrs->size() > 0 && (*fjPtrs)[0]->pt() > 400) {
      const FatJet& fj = *((*fjPtrs)[0]); 
      int inc1 = 0;
      if (DeltaR2(jets.cleaned_sorted[idx1]->base->eta(), 
                  jets.cleaned_sorted[idx1]->base->phi(),
                  fj.eta(), fj.phi()) < cfg.FATJETMATCHDR2) {
        inc1++;
      }
      if (DeltaR2(jets.cleaned_sorted[idx0]->base->eta(), 
                  jets.cleaned_sorted[idx0]->base->phi(),
                  fj.eta(), fj.phi()) < cfg.FATJETMATCHDR2) {
        if (inc1 == 0) {
          inc1 = 1;
          idx0 = idx1; 
        } else {
          idx0 = idx1 + inc1;
          inc1++;
        }
      }
      idx1 += inc1;
    } else if (gt.hbbm[shift] > 0) {
      if (gt.hbbjtidx[shift][0] == 0 || gt.hbbjtidx[shift][1] == 0) {
        idx0++; idx1++;
      }
      if (gt.hbbjtidx[shift][0] == 1 || gt.hbbjtidx[shift][1] == 1) {
        if (idx0 == 0) 
          idx1++;
        else {
          idx0++; idx1++;
        }
      }
    }
  }
  if (jets.cleaned.size() > idx1) {
    TLorentzVector v0 = jets.cleaned_sorted[idx0]->p4();
    TLorentzVector v1 = jets.cleaned_sorted[idx1]->p4();
    gt.jot12Mass[shift] = (v0 + v1).M();
    gt.jot12DPhi[shift] = v0.DeltaPhi(v1);
    gt.jot12DEta[shift] = fabs(v0.Eta() - v1.Eta());
  }
}

void HbbSystemOp::do_execute()
{
  auto& jets = **currentJES;
  int shift = jets.shift_idx;

  if (jets.central.size() < 2)
    return;

  btagsorted = jets.central; // copy
  sort(btagsorted.begin(), btagsorted.end(),
       analysis.useCMVA ?
        [](const JetWrapper *x, const JetWrapper *y) { return x->base->cmva > y->base->cmva; } :
        (analysis.useDeepCSV ?
         [](const JetWrapper *x, const JetWrapper *y) { return x->base->deepCSVb+x->base->deepCSVbb > y->base->deepCSVb+y->base->deepCSVbb; } :
         [](const JetWrapper *x, const JetWrapper *y) { return x->base->csv > y->base->csv; }
        )
      );
  map<const JetWrapper*, int> order; // needed for output indexing
  for (int i = 0; i != (int)jets.cleaned.size(); ++i)
    order[jets.cleaned[i]] = i;


  array<TLorentzVector,2> hbbd;
  for (int i = 0; i != 2; ++i)  {
    gt.hbbjtidx[shift][i] = order[btagsorted[i]];
    btagsorted[i]->p4(hbbd[i]);
  }

  TLorentzVector hbbsystem = hbbd[0] + hbbd[1];

  gt.hbbpt[shift] = hbbsystem.Pt();
  gt.hbbeta[shift] = hbbsystem.Eta();
  gt.hbbphi[shift] = hbbsystem.Phi();
  gt.hbbm[shift] = hbbsystem.M();

  array<TLorentzVector,2> hbbd_corr, hbbd_dcorr, hbbd_qcorr;
  if (gt.hbbm[shift] > 0) {
    for (int i = 0; i<2; i++) {
      int idx = gt.hbbjtidx[shift][i];
      (*hbbdJet) = jets.cleaned[idx];
      auto& hbbdJetRef = **hbbdJet;

      if (shift == jes2i(shiftjes::kNominal)) {
        deepreg->execute();
        gt.jotDeepBReg[i] = hbbdJetRef.breg;
        gt.jotDeepBRegWidth[i] = hbbdJetRef.bregwidth;
        gt.jotDeepBRegSampled[i] = (event.rng.normal() * hbbdJetRef.bregwidth)  + hbbdJetRef.breg;
      }
      auto scale_fn = [&](float x) { return 1 + (gt.jotPt[shift][idx] / hbbdJetRef.base->pt()) * (x - 1); };
      hbbd_dcorr[i].SetPtEtaPhiM(
            gt.jotPt[shift][idx] * scale_fn(gt.jotDeepBReg[i]),
            gt.jotEta[idx],
            gt.jotPhi[idx],
            gt.jotM[idx]
          );
      hbbd_qcorr[i].SetPtEtaPhiM(
            gt.jotPt[shift][idx] * scale_fn(gt.jotDeepBRegSampled[i]),
            gt.jotEta[idx],
            gt.jotPhi[idx],
            gt.jotM[idx]
          );

      // Shifted values for the jet energies to perform the b-jet regression
      if (shift == jes2i(shiftjes::kNominal)) {
        bdtreg->execute();
        gt.jotBReg[i] = hbbdJetRef.breg;
      }
      hbbd_corr[i].SetPtEtaPhiM(
            gt.jotBReg[i] * gt.jotPt[shift][idx],
            gt.jotEta[idx],
            gt.jotPhi[idx],
            gt.jotM[idx]
          );
    }

    TLorentzVector hbbsystem_corr = hbbd_corr[0] + hbbd_corr[1];
    gt.hbbm_reg[shift] = hbbsystem_corr.M();
    gt.hbbpt_reg[shift] = hbbsystem_corr.Pt();

    TLorentzVector hbbsystem_dcorr = hbbd_dcorr[0] + hbbd_dcorr[1];
    gt.hbbm_dreg[shift] = hbbsystem_dcorr.M();
    gt.hbbpt_dreg[shift] = hbbsystem_dcorr.Pt();

    TLorentzVector hbbsystem_qcorr = hbbd_qcorr[0] + hbbd_qcorr[1];
    gt.hbbm_qreg[shift] = hbbsystem_qcorr.M();
    gt.hbbpt_qreg[shift] = hbbsystem_qcorr.Pt();

  } // regression

  if (gt.hbbm > 0) {
    gt.hbbCosThetaJJ[shift] = hbbsystem.CosTheta();
    float csj1;
    if (analysis.bjetBDTReg) {
      if (hbbd_corr[0].Pt() > hbbd_corr[1].Pt())
        csj1 = CosThetaCollinsSoper(hbbd_corr[0], hbbd_corr[1]);
      else
        csj1 = CosThetaCollinsSoper(hbbd_corr[1], hbbd_corr[0]);
    } else {
      if (hbbd[0].Pt() > hbbd[1].Pt())
        csj1 = CosThetaCollinsSoper(hbbd[0], hbbd[1]);
      else
        csj1 = CosThetaCollinsSoper(hbbd[1], hbbd[0]);
    }
    gt.hbbCosThetaCSJ1[shift] = csj1;
  }

  if (gt.hbbm > 0 && gt.nLooseLep > 0 && shift <= jes2i(shiftjes::kJESTotalDown)) {
    TLorentzVector leptonP4, metP4, nuP4, *jet0P4{nullptr}, *jet1P4{nullptr}, WP4, topP4;
    float dRJet0W, dRJet1W;
    bool jet0IsCloser;
    leptonP4=(*looseLeps)[0]->p4();

    metP4.SetPtEtaPhiM(gt.pfmet[shift], 0, gt.pfmetphi[shift], 0);
    nuP4 = getNu4Momentum( leptonP4, metP4 );
    WP4 = leptonP4 + nuP4;

    // If using b-jet regression, use the regressed jets for the top mass reconstruction
    // Otherwise, use the un regressed jets
    if (analysis.bjetBDTReg) {
      jet0P4 = &hbbd_corr[0]; jet1P4 = &hbbd_corr[1];
    } else {
      jet0P4 = &hbbd[0]; jet1P4 = &hbbd[1];
    }

    dRJet0W=jet0P4->DeltaR(leptonP4);
    dRJet1W=jet1P4->DeltaR(leptonP4);
    jet0IsCloser = (dRJet0W < dRJet1W);

    topP4 = jet0IsCloser ? (*jet0P4)+WP4 : (*jet1P4)+WP4;
    gt.topMassLep1Met[shift] = topP4.M();
    gt.topWBosonCosThetaCS[shift] = CosThetaCollinsSoper(WP4, jet0IsCloser ? *jet0P4 : *jet1P4);
    gt.topWBosonPt  = WP4.Pt();
    gt.topWBosonEta = WP4.Eta();
    gt.topWBosonPhi = WP4.Phi();
  }
  if (gt.nLooseLep > 1 && gt.hbbm > 0) {
    TLorentzVector HP4;
    if (analysis.bjetBDTReg)
      HP4.SetPtEtaPhiM(gt.hbbpt_reg[0], gt.hbbeta[0], gt.hbbphi[0], gt.hbbm_reg[0]);
    else
      HP4.SetPtEtaPhiM(gt.hbbpt[0], gt.hbbeta[0], gt.hbbphi[0], gt.hbbm[0]);
    TLorentzVector ZHP4 = (*dilep) + HP4;
    gt.ZBosonLep1CosThetaStar = CosThetaStar(looseLeps->at(0)->p4(), looseLeps->at(1)->p4(), ZHP4);
  }
}
