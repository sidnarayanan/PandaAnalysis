#include "../interface/JetsMods.h"
#include "PandaAnalysis/Utilities/interface/Helicity.h"
#include "PandaAnalysis/Utilities/interface/NeutrinoSolver.h"
#include "TSystem.h"

using namespace pa;
using namespace std;
using namespace panda; 
using namespace fastjet;
using JECParams = JetCorrectorParameters;

inline float centralOnly(float x, float aeta, float def = -1) 
{
  return  aeta < 2.4 ? x : -1;
}

JetWrapper BaseJetMod::shiftJet(const Jet& jet, shiftjes shift, bool smear) 
{
  float pt = smear ? jet.ptSmear : jet.pt(); 
  if (shift != shiftjes::kNominal) { 
    int ishift = jes2i(shift);
    bool isUp = !(ishift % 2 == 0); 
    (*scaleUnc)[ishift]->setJetPt(pt);
    (*scaleUnc)[ishift]->setJetEta(jet.eta());
    pt *= (*scaleUnc)[ishift]->getUncertainty(isUp);
  }
  return JetWrapper(pt, jet);
}


void BaseJetMod::do_readData(TString dirPath) 
{
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
  scales["MC"] = new FactorizedJetCorrector(params);
  for (auto e : eraGroups) {
    basePath = dirPath+"/jec/"+jecVFull+"/"+campaign+"_"+jecReco+e+spacer+jecV;
    params = {
      JECParams((basePath+"_DATA_L1FastJet_"+jetType+".txt").Data()),
      JECParams((basePath+"_DATA_L2Relative_"+jetType+".txt").Data()),
      JECParams((basePath+"_DATA_L3Absolute_"+jetType+".txt").Data()),
      JECParams((basePath+"_DATA_L2L3Residual_"+jetType+".txt").Data())
    };
    scales["data"+e] = new FactorizedJetCorrector(params);
  }

  jer = new JERReader(dirPath+"/jec/"+jerV+"/"+jerV+"_MC_SF_"+jetType+".txt",
                      dirPath+"/jec/"+jerV+"/"+jerV+"_MC_PtResolution_"+jetType+".txt");


  basePath = dirPath+"/jec/"+jecVFull+"/"+campaign+"_";
  setScaleUnc("MC", (basePath+jecVFull+"_MC_UncertaintySources_"+jetType+".txt").Data());
  for (auto e : eraGroups) {
    setScaleUnc("data"+e, (basePath+jecReco+e+spacer+jecV+"_DATA_UncertaintySources_"+jetType+".txt").Data());
  }
}


void BaseJetMod::setScaleUnc(TString tag, TString path)
{
  scaleUncs[tag] = std::vector<JetCorrectionUncertainty*>(jes2i(shiftjes::N),nullptr);  
  JESLOOP {
    if (shift % 2 == 0)
      continue; 
    TString shiftName = jesName(i2jes(shift));
    shiftName.ReplaceAll("JES","");
    shiftName = shiftName(0,shiftName.Length()-2);
    scaleUncs[tag][shift] = new JetCorrectionUncertainty(JECParams(path.Data(), shiftName.Data()));
    scaleUncs[tag][shift+1] = scaleUncs[tag][shift]; // Down = Up 
  }
}

void JetMod::setupJES()
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
        scale = scales[iter.first];
        return;
      }
    }
  } else {
    scaleUnc = &(scaleUncs["MC"]);
    scale = scales["MC"];
  }
}


void JetMod::varyJES()
{
  JESLOOP {
    auto& jets = (*jesShifts)[shift];
    jets.reserve(ak4Jets->size());
    for (auto &j : *ak4Jets) {
      jets.all.push_back(shiftJet(j, i2jes(shift), analysis.hbb && !analysis.isData));
    }
    jets.sort();
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
    bool isNominal = (shift == jes2i(shiftjes::kNominal));
    bool metShift = (i2jes(shift) <= shiftjes::kJESTotalDown);
    JESHandler& jets = (*jesShifts)[shift];
    currentJES = &jets;
    for (auto& jw : jets.all) {
      currentJet = &jw;
      auto& jet = jw.get_base();
      float aeta = abs(jet.eta());
      float pt = jw.pt;
      if (aeta > maxJetEta || pt < minMinJetPt) 
        continue;
      if (isMatched(matchLeps,0.16,jet.eta(),jet.phi()))
        continue;
      if (!analysis.hbb && isMatched(matchPhos,0.16,jet.eta(),jet.phi()))
        continue;
      if ((analysis.vbf || analysis.hbb) && !jet.loose)
        continue;

      flavor->execute();

      float csv = centralOnly(analysis.year==2016 ? jet.csv : jet.deepCSVb, aeta);
      float cmva = centralOnly(jet.cmva, aeta);

      if (pt > cfg.minBJetPt && aeta < 2.4) { // b jets
        if (csvLoose(csv)) { 
          ++(gt.jetNBtags[shift]);
          if (csvMed(csv)) {
            ++(gt.jetNMBtags[shift]);
          }
        }

        jets.bcand.push_back(&jw); 
      } 


      if (pt > cfg.minJetPt) {
        // for H->bb, don't consider any jet based NJETSAVED, 
        // for other analyses, consider them, just don't save them
        if ((analysis.hbb || analysis.monoh) && (int)jets.cleaned.size() >= cfg.NJETSAVED)
          continue;

        jets.cleaned.push_back(&jw);
        
        if (jets.cleaned.size() < 3) {
          if (utils.getCorr(cBadECALJets, jet.eta(), jet.phi()) > 0)
            gt.badECALFilter = 0;
        }

        if (analysis.fatjet) {
          isojet->execute();
          if (jw.iso && csvLoose(csv))
            ++(gt.isojetNBtags[shift]);
        }

        if (aeta < 2.4) {
          jets.central.push_back(&jw);
          
          int njet = jets.central.size() - 1;
          gt.nJet[shift] = njet + 1;
          if (isNominal) {
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

        if (isNominal && aeta < 3.0) {
          gt.barrelHT += jet.pt();
          vBarrelJets += jet.p4();
        }

        int njet = jets.cleaned.size() - 1;
        if (njet < 2 || ((analysis.hbb || analysis.monoh) && njet < cfg.NJETSAVED)) {
          gt.jotPt[shift][njet] = pt;
          if (isNominal) {
            gt.jotEta[njet] = jet.eta();
            gt.jotPhi[njet] = jet.phi();
            gt.jotM[njet] = jet.m();
            gt.jotCSV[njet] = csv;
            gt.jotCMVA[njet] = cmva;
            gt.jotVBFID[njet] = (aeta < 2.4) ? (jet.monojet ? 1 : 0) : 1;

            bjetreg->execute();
          }
        }

        TLorentzVector vJet; vJet.SetPtEtaPhiM(pt, jet.eta(), jet.phi(), jet.m());
        if (metShift && njet < nJetDPhi) { // only do this for fully-correlated shifts
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
    if (metShift) 
      vbf->execute();
    hbb->execute();

  } // shift loop 

  gt.barrelHTMiss = vBarrelJets.Pt();
  
}


void JetFlavorMod::partonFlavor()
{
  auto& jw = **currentJet;
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


void JetFlavorMod::clusteredFlavor()
{
  auto& jw = **currentJet;
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


void JetFlavorMod::do_execute()
{
  if (analysis.jetFlavorPartons)
    partonFlavor();
  else
    clusteredFlavor();
}


void IsoJetMod::do_execute()
{
  auto& jw = **currentJet;
  auto& jets = **currentJES;
  auto& jet = jw.get_base();
  float maxIsoEta = analysis.monoh ? 4.5 : 2.5;
  bool isIsoJet = ( 
        gt.nFatjet == 0 || 
        (fabs(jet.eta()) < maxIsoEta && 
         DeltaR2(gt.fjEta,gt.fjPhi,jet.eta(),jet.phi()) > cfg.FATJETMATCHDR2) 
      ); 

  jw.iso = isIsoJet; 

  if (isIsoJet) 
    jets.iso.push_back(&jw);
}

void BJetRegMod::do_execute()
{
  auto& jw = **currentJet;
  auto& jet = jw.get_base();
  auto& jets = **currentJES; 

  int N = jets.cleaned.size() - 1;

  gt.jotEMF[N] = jet.cef + jet.nef;
  gt.jotHF[N] = jet.chf + jet.nhf;
  gt.jotLep1Pt[N] = 0;
  gt.jotTrk1Pt[N] = 0;
  gt.jotNLep[N] = 0;
  for (const panda::ConstRef<panda::PFCand> &c_iter : jet.constituents) {
    if (!c_iter.isValid())
      continue;
    auto *pf = c_iter.get();
    if (pf->q() != 0) {
      float pt = pf->pt();
      gt.jotTrk1Pt[N] = max(pt, gt.jotTrk1Pt[N]);
      int pdgid = abs(pf->pdgId());
      if (pdgid == 11 || pdgid == 13) {
        gt.jotNLep[N]++;
        if (pt > gt.jotLep1Pt[N]) {
          gt.jotLep1Pt[N] = pt;
          gt.jotLep1PtRel[N] = pf->p4().Perp(jet.p4().Vect());
          gt.jotLep1DeltaR[N] = sqrt(DeltaR2(pf->eta(), pf->phi(), jet.eta(), jet.phi()));
        }
      }
    }
  }

  auto& vert = jet.secondaryVertex;
  if (vert.isValid()) {
    gt.jotVtxPt[N] = vert->pt();
    gt.jotVtxMass[N] = vert->m();
    gt.jotVtx3DVal[N] = vert->vtx3DVal;
    gt.jotVtx3DErr[N] = vert->vtx3DeVal;
    gt.jotVtxNtrk[N] = vert->ntrk;
  }
}

void VBFSystemMod::do_execute()
{
  auto& jets = **currentJES; 

  int shift = jets.shift_idx;

  if (jets.cleaned.size() > 1) {
    TLorentzVector v0 = jets.cleaned_sorted[0]->p4();
    TLorentzVector v1 = jets.cleaned_sorted[1]->p4();
    gt.jot12Mass[shift] = (v0 + v1).M();
    gt.jot12DPhi[shift] = v0.DeltaPhi(v1);
    gt.jot12DEta[shift] = fabs(v0.Eta() - v1.Eta());
  }
}

void HbbSystemMod::do_readData(TString dirPath)
{
  delete bjetreg_vars; 
  delete bjetregReader; // we're going to recreate these guys 

  bjetreg_vars = new float[15];
  bjetregReader = new TMVA::Reader("!Color:!Silent");
  bjetregReader->AddVariable("jetPt[hbbjtidx[0]]",
                             &(bjetreg_vars[ 0]) );
  bjetregReader->AddVariable("jetEta[hbbjtidx[0]]",
                             &(bjetreg_vars[ 1]) );
  bjetregReader->AddVariable("jetLeadingTrkPt[hbbjtidx[0]]",
                              &(bjetreg_vars[ 2]) );
  bjetregReader->AddVariable("jetLeadingLepPt[hbbjtidx[0]]",
                             &(bjetreg_vars[ 3]) );
  bjetregReader->AddVariable("jetEMFrac[hbbjtidx[0]]",
                             &(bjetreg_vars[ 4]) );
  bjetregReader->AddVariable("jetHadFrac[hbbjtidx[0]]",
                             &(bjetreg_vars[ 5]) );
  bjetregReader->AddVariable("jetLeadingLepDeltaR[hbbjtidx[0]]",
                             &(bjetreg_vars[ 6]) );
  bjetregReader->AddVariable("jetLeadingLepPtRel[hbbjtidx[0]]",
                             &(bjetreg_vars[ 7]) );
  bjetregReader->AddVariable("jetvtxPt[hbbjtidx[0]]",
                             &(bjetreg_vars[ 8]) );
  bjetregReader->AddVariable("jetvtxMass[hbbjtidx[0]]",
                             &(bjetreg_vars[ 9]) );
  bjetregReader->AddVariable("jetvtx3Dval[hbbjtidx[0]]",
                             &(bjetreg_vars[10]) );
  bjetregReader->AddVariable("jetvtx3Derr[hbbjtidx[0]]",
                             &(bjetreg_vars[11]) );
  bjetregReader->AddVariable("jetvtxNtrk[hbbjtidx[0]]",
                             &(bjetreg_vars[12]) );
  bjetregReader->AddVariable("evalEt(jetPt[hbbjtidx[0]],jetEta[hbbjtidx[0]],jetPhi[hbbjtidx[0]],jetE[hbbjtidx[0]])",
                             &(bjetreg_vars[13]) );
  bjetregReader->AddVariable("evalMt(jetPt[hbbjtidx[0]],jetEta[hbbjtidx[0]],jetPhi[hbbjtidx[0]],jetE[hbbjtidx[0]])",
                             &(bjetreg_vars[14]) );

  gSystem->Exec(
      Form("wget -nv -O %s/trainings/bjetregression.weights.xml http://t3serv001.mit.edu/~dhsu/pandadata/trainings/bjet_regression_v1_fromBenedikt.weights.xml",dirPath.Data())
    );
  bjetregReader->BookMVA("BDT method", 
                         dirPath+"trainings/bjetregression.weights.xml" );

}

void HbbSystemMod::do_execute()
{
  auto& jets = **currentJES; 
  int shift = jets.shift_idx;

  if (jets.central.size() < 2)
    return;

  btagsorted = jets.central; // copy
  sort(btagsorted.begin(), btagsorted.end(),
       analysis.useCMVA ? 
        [](const JetWrapper *x, const JetWrapper *y) { return x->base->cmva > y->base->cmva; } :   
        [](const JetWrapper *x, const JetWrapper *y) { return x->base->csv > y->base->csv; }
      );
  map<const JetWrapper*, int> order; // needed for output indexing
  for (int i = 0; i != (int)jets.cleaned.size(); ++i) 
    order[jets.cleaned[i]] = i;
 

  gt.hbbjtidx[shift][0] = order[btagsorted[0]];
  gt.hbbjtidx[shift][1] = order[btagsorted[1]];

  vector<TLorentzVector> hbbd(2);
  for (int i = 0; i != 2; ++i) {
    btagsorted[i]->p4(hbbd[i]);
  }
  TLorentzVector hbbsystem = hbbd[0] + hbbd[1];

  gt.hbbpt[shift] = hbbsystem.Pt();
  gt.hbbeta[shift] = hbbsystem.Eta();
  gt.hbbphi[shift] = hbbsystem.Phi();
  gt.hbbm[shift] = hbbsystem.M();

  array<TLorentzVector,2> hbbd_corr;
  if (analysis.bjetRegression && gt.hbbm[shift] > 0) {
    for (int i = 0; i<2; i++) {
      int idx = gt.hbbjtidx[shift][i];
      // Shifted values for the jet energies to perform the b-jet regression
      bjetreg_vars[0] = gt.jotPt[shift][idx];
      bjetreg_vars[1] = gt.jotEta[idx];
      bjetreg_vars[2] = gt.jotTrk1Pt[idx];
      bjetreg_vars[3] = gt.jotLep1Pt[idx];
      bjetreg_vars[4] = gt.jotEMF[idx];
      bjetreg_vars[5] = gt.jotHF[idx];
      bjetreg_vars[6] = gt.jotLep1DeltaR[idx];
      bjetreg_vars[7] = gt.jotLep1PtRel[idx];
      bjetreg_vars[8] = gt.jotVtxPt[idx];
      bjetreg_vars[9] = gt.jotVtxMass[idx];
      bjetreg_vars[10]= gt.jotVtx3DVal[idx];
      bjetreg_vars[11]= gt.jotVtx3DErr[idx];
      bjetreg_vars[12]= gt.jotVtxNtrk[idx];
      bjetreg_vars[13]= hbbd[i].Et();
      bjetreg_vars[14]= hbbd[i].Mt();

      gt.jotBReg[shift][i] = (bjetregReader->EvaluateRegression("BDT method"))[0];
      hbbd_corr[i].SetPtEtaPhiM(
            gt.jotBReg[shift][i] * gt.jotPt[shift][idx],
            gt.jotEta[idx],
            gt.jotPhi[idx],
            gt.jotM[idx]
          );
    }

    TLorentzVector hbbsystem_corr = hbbd_corr[0] + hbbd_corr[1];
    gt.hbbm_reg[shift] = hbbsystem_corr.M(); 
    gt.hbbpt_reg[shift] = hbbsystem_corr.Pt();

  } // regression

  if (gt.hbbm > 0) {
    gt.hbbCosThetaJJ[shift] = hbbsystem.CosTheta();
    float csj1;
    if (analysis.bjetRegression) {
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
    if (analysis.bjetRegression) {
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
    if (analysis.bjetRegression)
      HP4.SetPtEtaPhiM(gt.hbbpt_reg[0], gt.hbbeta[0], gt.hbbphi[0], gt.hbbm_reg[0]);
    else
      HP4.SetPtEtaPhiM(gt.hbbpt[0], gt.hbbeta[0], gt.hbbphi[0], gt.hbbm[0]);
    TLorentzVector ZHP4 = (*dilep) + HP4;
    gt.ZBosonLep1CosThetaStar = CosThetaStar(looseLeps->at(0)->p4(), looseLeps->at(1)->p4(), ZHP4);
  }
}

