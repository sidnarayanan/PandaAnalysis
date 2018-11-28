#include "../interface/FatJetsOps.h"
#include "PandaAnalysis/Utilities/interface/Helicity.h"
#include "PandaAnalysis/Utilities/interface/EnergyCorrelations.h"
#include "fastjet/contrib/Njettiness.hh"

using namespace pa;
using namespace std;
using namespace panda;
using namespace fastjet;

void FatJetOp::setupJES()
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
        return;
      }
    }
  } else {
    scaleUnc = &(scaleUncs["MC"]);
  }
}

void FatJetOp::do_execute()
{
  setupJES();

  gt.nFatJet=0;
  int fatjet_counter=-1;
  for (auto* fjPtr : *fjPtrs) {
    auto& fj = *fjPtr;
    ++fatjet_counter;
    bool doSmear = !analysis.isData && recalcJER;
    JetWrapper jwNominal = shiftJet(fj, shiftjes::kNominal, doSmear); 
    float pt = jwNominal.pt;
    float rawpt = fj.rawPt;
    float eta = fj.eta();
    float mass = fj.m();
    float ptcut = 200;
    if (analysis.deep)
      ptcut = 400;

    float bestPt = pt;
    if (analysis.varyJES || analysis.varyJESTotal) {
      bestPt = TMath::Max(bestPt, fj.ptCorrUp);
      bestPt = TMath::Max(bestPt, fj.ptCorrDown);
    }
    if (bestPt<ptcut || fabs(eta)>2.4 || !fj.monojet)
      continue;

    float phi = fj.phi();
    if (isMatched(matchLeps.get(),cfg.FATJETMATCHDR2,eta,phi) ||
        isMatched(matchPhos.get(),cfg.FATJETMATCHDR2,eta,phi)) {
      continue;
    }

    int iFJ = gt.nFatJet;
    double best_doubleB = 0; 
    gt.nFatJet++;
    if (iFJ < nMaxFJ) {
      gt.nFatJetTrunc++; 
      gt.fjIsClean[iFJ] = fatjet_counter==0 ? 1 : 0;
      gt.fjEta[iFJ] = eta;
      gt.fjPhi[iFJ] = phi;
      gt.fjRawPt[iFJ] = rawpt;
      float corrweight = getMSDCorr(pt,eta);
      JESLOOP {
        JetWrapper jw = shift == 0 ? jwNominal : shiftJet(fj, i2jes(shift), doSmear);
        gt.fjPt[shift][iFJ] = jw.pt;
        gt.fjM[shift][iFJ] = mass * jw.scale();
        gt.fjMSD[shift][iFJ] = fj.mSD * jw.scale();
        gt.fjMSD_corr[shift][iFJ] = corrweight*gt.fjMSD[shift][iFJ];
      }

      if (analysis.recalcECF && event.recoil.max < 175) {
        substructure->run(fj);
      }

      // now we do substructure
      gt.fjTau32[iFJ] = clean(fj.tau3/fj.tau2);
      gt.fjTau32SD[iFJ] = clean(fj.tau3SD/fj.tau2SD);
      gt.fjTau21[iFJ] = clean(fj.tau2/fj.tau1);
      gt.fjTau21SD[iFJ] = clean(fj.tau2SD/fj.tau1SD);

      for (auto ibeta : cfg.ibetas) {
        for (auto N : cfg.Ns) {
          for (auto order : cfg.orders) {
            GeneralTree::ECFParams p;
            p.order = order; p.N = N; p.ibeta = ibeta; 
            if (gt.fjIsClean[iFJ] || true)
              gt.fjECFNs[p][iFJ] = fj.get_ecf(order,N,ibeta);
            else
              gt.fjECFNs[p][iFJ] = fj.get_ecf(order,N,ibeta);
          }
        }
      } //loop over betas
      gt.fjHTTMass[iFJ] = fj.htt_mass;
      gt.fjHTTFRec[iFJ] = fj.htt_frec;

      if (analysis.vqqhbb) {
        // jet charge
        gt.fjQ[iFJ] = 0;
        for (const auto& pf : fj.constituents) {
          gt.fjQ[iFJ] += pf->q() * pf->pt() / jwNominal.pt;
        }
      }

      std::vector<MicroJet const*> subjets;
      for (int iS(0); iS != fj.subjets.size(); ++iS)
        subjets.push_back(&fj.subjets.objAt(iS));


      auto& mya = analysis; // local scope
      auto mycsv = [&mya](const MicroJet& j) { return (mya.year == 2016 ? j.csv : j.deepCSVb+j.deepCSVbb); };
      auto csvsort = [&mycsv](MicroJet const* j1, MicroJet const* j2) -> bool {
              return mycsv(*j1) > mycsv(*j2);
            };

      std::sort(subjets.begin(),subjets.end(),csvsort);
      if (subjets.size() > 0) {
        gt.fjMaxCSV[iFJ] = mycsv(*(subjets[0]));
        gt.fjMinCSV[iFJ] = mycsv(*(subjets.back()));
        if (subjets.size()>1) {
          gt.fjSubMaxCSV[iFJ] = mycsv(*(subjets[1]));
        }
      }

      gt.fjDoubleCSV[iFJ] =  fj.double_sub;
      gt.fjDeepProbH[iFJ] = fj.deepBBprobH;
      gt.fjDeepProbbb[iFJ] = fj.deepCSVbb; 
      if (gt.fjDoubleCSV[iFJ] > best_doubleB) {
        best_doubleB = gt.fjDoubleCSV[iFJ]; 
        gt.fjHiggsIdx = iFJ; 
      }

      if (analysis.monoh && iFJ == 0) {
        for (int iSJ=0; iSJ!=fj.subjets.size(); ++iSJ) {
          auto& subjet = fj.subjets.objAt(iSJ);
          gt.fjsjPt[iSJ]=subjet.pt();
          gt.fjsjEta[iSJ]=subjet.eta();
          gt.fjsjPhi[iSJ]=subjet.phi();
          gt.fjsjM[iSJ]=subjet.m();
          gt.fjsjCSV[iSJ]=mycsv(subjet);
          gt.fjsjQGL[iSJ]=subjet.qgl;
        }
      }
      if (!analysis.isData && fj.matchedGenJet.isValid())
        gt.fjGenNumB[iFJ] = fj.matchedGenJet.get()->numB;
      else
        gt.fjGenNumB[iFJ] = 0;
    }
  }
  gt.fjVIdx = 1 - gt.fjHiggsIdx; 

  if (analysis.hbb) {
    if (gt.nFatJet > 0 && gt.nLooseLep > 1) {
      TLorentzVector HP4;
      HP4.SetPtEtaPhiM(gt.fjPt[0][0], gt.fjEta[0], gt.fjPhi[0], gt.fjMSD_corr[0][0]);
      TLorentzVector ZHP4 = (*dilep) + HP4;
      gt.ZBosonLep1CosThetaStarFJ = CosThetaStar((*looseLeps)[0]->p4(), (*looseLeps)[1]->p4(), ZHP4);
    }
  }
}


float FatJetOp::getMSDCorr(float puppipt, float puppieta)
{

  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;

  genCorr = utils.puppisd_corrGEN->Eval( puppipt );
  if (fabs(puppieta) <= 1.3) {
    recoCorr = utils.puppisd_corrRECO_cen->Eval(puppipt);
  } else {
    recoCorr = utils.puppisd_corrRECO_for->Eval(puppipt);
  }
  //printf("genCorr=%.3f, recoCorr=%.3f\n",genCorr,recoCorr);
  totalWeight = genCorr * recoCorr;

  return totalWeight;
}

const GenParticle* FatJetMatchingOp::matchGen(double eta, double phi, double radius, int pdgid) const
{
  const GenParticle* found=NULL;
  double r2 = radius*radius;
  pdgid = abs(pdgid);

  for (auto iG=genObjects.begin();
      iG!=genObjects.end(); ++iG) {
    if (found!=NULL)
      break;
    if (pdgid!=0 && abs(iG->first->pdgid)!=pdgid)
      continue;
    if (DeltaR2(eta,phi,iG->first->eta(),iG->first->phi())<r2)
      found = iG->first;
  }

  return found;
}


void FatJetMatchingOp::do_execute()
{
  if (fjPtrs->size() == 0)
    return; 

  int pdgidTarget=0;
  if (!analysis.isData && analysis.processType>=kTT && analysis.processType<=kSignal) {
    switch(analysis.processType) {
      case kTop:
      case kTT:
      case kSignal:
        pdgidTarget=6;
        break;
      case kV:
        pdgidTarget=24;
        break;
      case kH:
        pdgidTarget=25;
        break;
      default:
        // analysis.processType>=kTT means we should never get here
        logger.error("FatJetMatchingOp::do_execute","Reached an unknown process type");
    }

    std::vector<int> targets;

    int nGen = genP->size();
    for (int iG=0; iG!=nGen; ++iG) {
      auto& part = pToGRef((*genP)[iG]);
      int pdgid = part.pdgid;
      int abspdgid = abs(pdgid);
      if (abspdgid == pdgidTarget)
        targets.push_back(iG);
    } //looking for targets

    for (int iG : targets) {
      auto& part = pToGRef((*genP)[iG]);

      // check there is no further copy:
      bool isLastCopy=true;
      for (int jG : targets) {
        if (pToGPtr((*genP)[jG])->parent.get() == &part) {
          isLastCopy=false;
          break;
        }
      }
      if (!isLastCopy)
        continue;

      // (a) check it is a hadronic decay and if so, (b) calculate the size
      if (analysis.processType==kTop||analysis.processType==kTT) {

        // first look for a W whose parent is the top at iG, or a W further down the chain
        const GenParticle* lastW(0);
        for (int jG=0; jG!=nGen; ++jG) {
          const GenParticle& partW = pToGRef((*genP)[jG]);
          if (TMath::Abs(partW.pdgid)==24 && partW.pdgid*part.pdgid>0) {
            // it's a W and has the same sign as the top
            if (!lastW && partW.parent.get() == &part) {
              lastW = &partW;
            } else if (lastW && partW.parent.get() == lastW) {
              lastW = &partW;
            }
          }
        } // looking for W
        if (!lastW) {// ???
          continue;
        }
        auto& partW(*lastW);

        // now look for b or W->qq
        int iB=-1, iQ1=-1, iQ2=-1;
        double size=0, sizeW=0;
        for (int jG=0; jG!=nGen; ++jG) {
          auto& partQ = pToGRef((*genP)[jG]);
          int pdgidQ = partQ.pdgid;
          int abspdgidQ = TMath::Abs(pdgidQ);
          if (abspdgidQ>5)
            continue;
          if (abspdgidQ==5 && iB<0 && partQ.parent.get() == &part) {
            // only keep first copy
            iB = jG;
            size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),size);
          } else if (abspdgidQ<5 && partQ.parent.get() == &partW) {
            if (iQ1<0) {
              iQ1 = jG;
              size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),
                  size);
              sizeW = TMath::Max(DeltaR2(partW.eta(),partW.phi(),partQ.eta(),partQ.phi()),
                  sizeW);
            } else if (iQ2<0) {
              iQ2 = jG;
              size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),
                  size);
              sizeW = TMath::Max(DeltaR2(partW.eta(),partW.phi(),partQ.eta(),partQ.phi()),
                  sizeW);
            }
          }
          if (iB>=0 && iQ1>=0 && iQ2>=0)
            break;
        } // looking for quarks


        bool isHadronic = (iB>=0 && iQ1>=0 && iQ2>=0); // all 3 quarks were found
        if (isHadronic)
          genObjects[&part] = size;

        bool isHadronicW = (iQ1>=0 && iQ2>=0);
        if (isHadronicW)
          genObjects[&partW] = sizeW;

      } else { // these are W,Z,H - 2 prong decays

        int iQ1=-1, iQ2=-1;
        double size=0;
        for (int jG=0; jG!=nGen; ++jG) {
          auto& partQ = pToGRef((*genP)[jG]);
          int pdgidQ = partQ.pdgid;
          int abspdgidQ = TMath::Abs(pdgidQ);
          if (abspdgidQ>5)
            continue;
          if (partQ.parent.get() == &part) {
            if (iQ1<0) {
              iQ1=jG;
              size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),
                  size);
            } else if (iQ2<0) {
              iQ2=jG;
              size = TMath::Max(DeltaR2(part.eta(),part.phi(),partQ.eta(),partQ.phi()),
                  size);
            }
          }
          if (iQ1>=0 && iQ2>=0)
            break;
        } // looking for quarks

        bool isHadronic = (iQ1>=0 && iQ2>=0); // both quarks were found

        // add to collection
        if (isHadronic)
          genObjects[&part] = size;
      }

    } // loop over targets
  } // process is interesting

  int iFJ = -1; 
  for (auto* fj : *fjPtrs) {
    ++iFJ; 
    // first see if jet is matched
    auto* matched = matchGen(fj->eta(),fj->phi(),1.5,pdgidTarget);
    if (matched!=nullptr) {
      gt.fjIsMatched[iFJ] = 1;
      gt.fjGenPt[iFJ] = matched->pt();
      gt.fjGenSize[iFJ] = genObjects[matched];
    } else {
      gt.fjIsMatched[iFJ] = 0;
    }
    if (pdgidTarget==6) { // matched to top; try for W
      auto* matchedW = matchGen(fj->eta(),fj->phi(),1.5,24);
      if (matchedW!=nullptr) {
        gt.fjIsWMatched[iFJ] = 1;
        gt.fjGenWPt[iFJ] = matchedW->pt();
        gt.fjGenWSize[iFJ] = genObjects[matchedW];
      } else {
        gt.fjIsWMatched[iFJ] = 0;
      }
    }

    bool found_b_from_g=false;
    int bs_inside_cone=0;
    int has_gluon_splitting=0;
    const GenParticle* first_b_mo(0);
    // now get the highest pT gen particle inside the jet cone
    for (auto* genptr : *genP) {
      auto& gen = pToGRef(genptr);
      float pt = gen.pt();
      int pdgid = gen.pdgid;
      if (pt>(gt.fjHighestPtGenPt[iFJ])
          && DeltaR2(gen.eta(),gen.phi(),fj->eta(),fj->phi())<cfg.FATJETMATCHDR2) {
        gt.fjHighestPtGenPt[iFJ] = pt;
        gt.fjHighestPtGen[iFJ] = pdgid;
      }

      if (gen.parent.isValid() && gen.parent->pdgid==gen.pdgid)
        continue;

      //count bs and cs
      int apdgid = abs(pdgid);
      if (apdgid!=5 && apdgid!=4)
        continue;

      if (DeltaR2(gen.eta(),gen.phi(),fj->eta(),fj->phi())<cfg.FATJETMATCHDR2) {
        gt.fjNHF[iFJ]++;
        if (apdgid==5) {
          if (gen.parent.isValid() && gen.parent->pdgid==21 && gen.parent->pt()>20) {
            if (!found_b_from_g) {
              found_b_from_g=true;
              first_b_mo=gen.parent.get();
              bs_inside_cone+=1;
            } else if (gen.parent.get()==first_b_mo) {
              bs_inside_cone+=1;
              has_gluon_splitting=1;
            } else {
              bs_inside_cone+=1;
            }
          } else {
            bs_inside_cone+=1;
          }
        }
      }
    }

    gt.fjNbs[iFJ]=bs_inside_cone;
    gt.fjgbb[iFJ]=has_gluon_splitting;

    if (analysis.btagSFs && iFJ == 0) {
      // now get the subjet btag SFs
      vector<btagcand> sj_btagcands;
      vector<double> sj_sf_cent, sj_sf_bUp, sj_sf_bDown, sj_sf_mUp, sj_sf_mDown;
      int nSJ = fj->subjets.size();
      for (int iSJ=0; iSJ!=nSJ; ++iSJ) {
        auto& subjet = fj->subjets.objAt(iSJ);
        int flavor=0;
        for (auto* genptr : *genP) {
          auto& gen = pToGRef(genptr);
          int apdgid = abs(gen.pdgid);
          if (apdgid==0 || (apdgid>5 && apdgid!=21)) // light quark or gluon
            continue;
          double dr2 = DeltaR2(subjet.eta(),subjet.phi(),gen.eta(),gen.phi());
          if (dr2<0.09) {
            if (apdgid==4 || apdgid==5) {
              flavor=apdgid;
              break;
            } else {
              flavor=0;
            }
          }
        } // finding the subjet flavor

        float pt = subjet.pt();
        float btagUncFactor = 1;
        float eta = subjet.eta();
        double eff(1),sf(1),sfUp(1),sfDown(1);
        if (flavor==5) {
          eff = utils.getCorr(cCSVBL, pt, fabs(eta));
        } else if (flavor==4) {
          eff = utils.getCorr(cCSVCL, pt, fabs(eta));
        } else {
          eff = utils.getCorr(cCSVLL, pt, fabs(eta));
        }
        if (analysis.hbb)
          utils.btag->calcSF(bSubJetM,flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
        else
          utils.btag->calcSF(bSubJetL,flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
        sj_btagcands.push_back(btagcand(iSJ,flavor,eff,sf,sfUp,sfDown));
        sj_sf_cent.push_back(sf);
        if (flavor>0) {
          sj_sf_bUp.push_back(sfUp); sj_sf_bDown.push_back(sfDown);
          sj_sf_mUp.push_back(sf); sj_sf_mDown.push_back(sf);
        } else {
          sj_sf_bUp.push_back(sf); sj_sf_bDown.push_back(sf);
          sj_sf_mUp.push_back(sfUp); sj_sf_mDown.push_back(sfDown);
        }

      } // loop over subjets
      utils.btag->evalSF(sj_btagcands,sj_sf_cent,GeneralTree::bCent,GeneralTree::bSubJet,true);
      utils.btag->evalSF(sj_btagcands,sj_sf_bUp,GeneralTree::bBUp,GeneralTree::bSubJet,true);
      utils.btag->evalSF(sj_btagcands,sj_sf_bDown,GeneralTree::bBDown,GeneralTree::bSubJet,true);
      utils.btag->evalSF(sj_btagcands,sj_sf_mUp,GeneralTree::bMUp,GeneralTree::bSubJet,true);
      utils.btag->evalSF(sj_btagcands,sj_sf_mDown,GeneralTree::bMDown,GeneralTree::bSubJet,true);
    }

  }
}

float HRTagOp::getMSDCorr(float puppipt, float puppieta)
{

  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;

  genCorr = utils.puppisd_corrGEN->Eval( puppipt );
  if (fabs(puppieta) <= 1.3) {
    recoCorr = utils.puppisd_corrRECO_cen->Eval(puppipt);
  } else {
    recoCorr = utils.puppisd_corrRECO_for->Eval(puppipt);
  }
  totalWeight = genCorr * recoCorr;

  return totalWeight;
}


void HRTagOp::do_execute()
{
  // first loop through genP and find any partons worth saving
  // if this is a signal
  int i_parton = 0; 
  float minPt = 0;
  if (analysis.processType == kTop || analysis.processType == kTT) {
    for (auto *pptr : *genP) {
      auto& p = pToGRef(pptr);
      if (abs(p.pdgid) != 6)
        continue;
      float pt = p.pt();
      if (pt < minPt)
        continue;
      float eta = p.eta(), phi = p.phi();
      if (hasChild(p, *genP))
        continue;
      const GenParticle *W{nullptr}, *b{nullptr}, *q0{nullptr}, *q1{nullptr};
      // first loop through: find W and b
      for (auto *childptr : *genP) {
        auto& child = pToGRef(childptr);
        int id = abs(child.pdgid);
        if (id != 5 && id != 24)
          continue;
        if (hasChild(child, *genP))
          continue;
        if (!isAncestor(child, p))
          continue;
        if (id == 5) {
          b = &child;
        } else {
          W = &child;
        }
      } // W and b loop
      if (W == nullptr || b == nullptr)
        continue;
      // second loop through: find qq'
      for (auto *childptr : *genP) {
        auto& child = pToGRef(childptr);
        if (!child.parent.isValid() || child.parent.get() != W)
          continue;
        int id = abs(child.pdgid);
        if (id > 5)
          continue;
        if (q0 == nullptr)
          q0 = &child;
        else
          q1 = &child;
        if (q0 != nullptr && q1 != nullptr)
          break;
      } // qq' loop
      if (q0 == nullptr || q1 == nullptr)
        continue;
      // now fill the tree
      gt.i_evt = event.eventNumber;
      gt.i_parton = i_parton++;
      gt.npv = event.npv;
      gt.sampleType = 2;
      gt.gen_pt = pt; gt.gen_eta = eta; gt.gen_phi = p.phi();
      gt.gen_pdgid = p.pdgid;
      gt.gen_size = max({DeltaR2(eta, phi, b->eta(), b->phi()),
                         DeltaR2(eta, phi, q0->eta(), q0->phi()),
                         DeltaR2(eta, phi, q1->eta(), q1->phi())});
      // now find a good fat jet
      for (auto* fj : *fjPtrs) {
        if (DeltaR2(eta, phi, fj->eta(), fj->phi()) > 0.36)
          continue;
        fillJet(*fj);
        break;
      }
      gt.Fill();
      gt.Reset();
    }
  } else {
    for (auto *pptr : *genP) {
      auto& p = pToGRef(pptr);
      if (abs(p.pdgid) > 5 && abs(p.pdgid) != 21)
        continue;
      float pt = p.pt();
      if (pt < minPt)
        continue;
      float eta = p.eta(), phi = p.phi();
      if (!hard(p) || hasChild(p, *genP, true))
        continue;
      // now fill the tree
      gt.i_evt = event.eventNumber;
      gt.i_parton = i_parton++;
      gt.npv = event.npv;
      gt.sampleType = 0;
      gt.gen_pt = pt; gt.gen_eta = eta; gt.gen_phi = p.phi();
      gt.gen_pdgid = p.pdgid;
      gt.gen_size = 0;
      // now find a good fat jet
      for (auto* fj : *fjPtrs) {
        if (DeltaR2(eta, phi, fj->eta(), fj->phi()) > 0.36)
          continue;
        fillJet(*fj);
        break;
      }
      gt.Fill();
      gt.Reset();
    }
  }
}

void SubRunner::run(panda::FatJet& fj)
{
  VPseudoJet particles = convertPFCands(fj.constituents,doPuppi,0.01);

  ClusterSequenceArea seq(particles,*jetDef,*(utils.areaDef));
  VPseudoJet allJets(seq.inclusive_jets(0.));
  PseudoJet *pj{nullptr};
  double minDR2 = 999;
  for (auto &jet : allJets) {
    double dr2 = DeltaR2(jet.eta(),jet.phi_std(),fj.eta(),fj.phi());
    if (dr2<minDR2) {
      minDR2 = dr2;
      pj = &jet;
    }
  }
  if (pj == nullptr)
    return;
  VPseudoJet constituents = sorted_by_pt(pj->constituents());

  PseudoJet sd = (*utils.softDrop)(*pj);
  VPseudoJet sdConstituents = sorted_by_pt(sd.constituents());
  fj.mSD = sd.m();

  fj.tau1 = tauN->getTau(1, constituents);
  fj.tau2 = tauN->getTau(2, constituents);
  fj.tau3 = tauN->getTau(3, constituents);
  fj.tau1SD = tauN->getTau(1, sdConstituents);
  fj.tau2SD = tauN->getTau(2, sdConstituents);
  fj.tau3SD = tauN->getTau(3, sdConstituents);

  // get ecfs
  if (doECF) {
    unsigned nFilter = min(100, (int)sdConstituents.size());
    VPseudoJet sdConstsFiltered(sdConstituents.begin(), sdConstituents.begin() + nFilter);

    ecfcalc->calculate(sdConstsFiltered);
    for (auto iter = ecfcalc->begin(); iter != ecfcalc->end(); ++iter) {
      int N = iter.get<pa::ECFCalculator::nP>();
      int o = iter.get<pa::ECFCalculator::oP>();
      int beta = iter.get<pa::ECFCalculator::bP>();
      float ecf(iter.get<pa::ECFCalculator::ecfP>());
      if (!fj.set_ecf(o+1,N+1,beta,ecf)) {
        logger.error("SubRunner::run", Form("Failed to set ecf at %i %i %i\n", o+1, N+1, beta));
        throw std::runtime_error("");
      }
    }
  }

  // HTT
  PseudoJet httJet = htt->result(*pj);
  if (httJet != 0) {
    auto* s(static_cast<httwrapper::HEPTopTaggerV2Structure*>(httJet.structure_non_const_ptr()));
    fj.htt_mass = s->top_mass();
    fj.htt_frec = s->fRec();
  } else {
    fj.htt_mass = 0;
    fj.htt_frec = 0;
  }
}

void HRTagOp::fillJet(panda::FatJet& fj)
{
  if (analysis.reclusterFJ) 
    substructure->run(fj);

  gt.recoil = event.recoil.max;
  gt.clf_IsMatched = 1;
  gt.clf_Pt = fj.pt();
  gt.clf_Eta = fj.eta();
  gt.clf_Phi = fj.phi();
  gt.clf_M = fj.m();
  gt.clf_Tau32 = clean(fj.tau3 / fj.tau2);
  gt.clf_Tau21 = clean(fj.tau2 / fj.tau1);
  gt.clf_Tau32SD = clean(fj.tau3SD / fj.tau2SD);
  gt.clf_Tau21SD = clean(fj.tau2SD / fj.tau1SD);
  gt.clf_HTTFRec = fj.htt_frec;
  gt.clf_MSD = fj.mSD;
  gt.clf_MSD_corr = fj.mSD * getMSDCorr(fj.pt(),fj.eta());

  std::vector<MicroJet const*> subjets;
  for (int iS(0); iS != fj.subjets.size(); ++iS)
    subjets.push_back(&fj.subjets.objAt(iS));
  auto csvsort = [](MicroJet const* j1, MicroJet const* j2) -> bool {
          return j1->csv > j2->csv;
        };
  std::sort(subjets.begin(),subjets.end(),csvsort);
  if (subjets.size()>0)
    gt.clf_MaxCSV = subjets[0]->csv;

  for (auto ibeta : gt.get_ibetas()) {
    for (auto N : gt.get_Ns()) {
      for (auto order : gt.get_orders()) {
        HeavyResTree::ECFParams p;
        p.order = order; p.N = N; p.ibeta = ibeta;
        gt.clf_ECFNs[p] = fj.get_ecf(order,N,ibeta);
      }
    }
  } //loop ov
}
