#include "../interface/FatJetsMods.h"
#include "PandaAnalysis/Utilities/interface/Helicity.h"

using namespace pa;
using namespace std;
using namespace panda; 
namespace fj = fastjet;

void FatJetMod::setupJES()
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

void FatJetMod::do_execute()
{
  setupJES();

  gt.nFatjet=0;
  int fatjet_counter=-1;
  for (auto& fj : fatjets) {
    ++fatjet_counter;
    float pt = (analysis.hbb && !analysis.isData) ? fj.pt() : fj.ptSmear;
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
    if (isMatched(matchLeps.get(),cfg.FATJETMATCHDR2,eta,phi) || 
        isMatched(matchPhos.get(),cfg.FATJETMATCHDR2,eta,phi)) {
      continue;
    }

    gt.nFatjet++;
    if (gt.nFatjet==1) {
      *fj1 = &fj;
      gt.fjIsClean = fatjet_counter==0 ? 1 : 0;
      gt.fjEta = eta;
      gt.fjPhi = phi;
      gt.fjRawPt = rawpt;
      JESLOOP {
        JetWrapper jw = shiftJet(fj, i2jes(shift), analysis.hbb && !analysis.isData);
        gt.fjPt[shift] = jw.pt;
        gt.fjM[shift] = mass * jw.scale();
        gt.fjMSD[shift] = fj.mSD * jw.scale();
      }

      // mSD correction
      float corrweight = getMSDCorr(pt,eta);
      JESLOOP { 
        gt.fjMSD_corr[shift] = corrweight*gt.fjMSD[shift];
      }

      // now we do substructure
      gt.fjTau32 = clean(fj.tau3/fj.tau2);
      gt.fjTau32SD = clean(fj.tau3SD/fj.tau2SD);
      gt.fjTau21 = clean(fj.tau2/fj.tau1);
      gt.fjTau21SD = clean(fj.tau2SD/fj.tau1SD);

      for (auto ibeta : cfg.ibetas) {
        for (auto N : cfg.Ns) {
          for (auto order : cfg.orders) {
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

      std::vector<MicroJet const*> subjets;
      for (int iS(0); iS != fj.subjets.size(); ++iS)
        subjets.push_back(&fj.subjets.objAt(iS));


      auto& mya = analysis; // local scope
      auto mycsv = [&mya](const MicroJet& j) { return (mya.year == 2016 ? j.csv : j.deepCSVb); };
      auto csvsort = [&mycsv](MicroJet const* j1, MicroJet const* j2) -> bool {
              return mycsv(*j1) > mycsv(*j2);
            };

      std::sort(subjets.begin(),subjets.end(),csvsort);
      if (subjets.size()>0) {
        gt.fjMaxCSV = mycsv(*(subjets[0]));
        gt.fjMinCSV = mycsv(*(subjets.back()));
        if (subjets.size()>1) {
          gt.fjSubMaxCSV = mycsv(*(subjets[1]));
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
          gt.fjsjCSV[iSJ]=mycsv(subjet);
          gt.fjsjQGL[iSJ]=subjet.qgl;
        }
      }
    }
    if (!analysis.isData && fj.matchedGenJet.isValid())
      gt.fjGenNumB = fj.matchedGenJet.get()->numB;
    else 
      gt.fjGenNumB = 0;
  }

  if (analysis.hbb) {
    if (gt.nFatjet > 0 && gt.nLooseLep > 1) {
      TLorentzVector HP4;
      HP4.SetPtEtaPhiM(gt.fjPt[0], gt.fjEta, gt.fjPhi, gt.fjMSD_corr[0]);
      TLorentzVector ZHP4 = (*dilep) + HP4;
      gt.ZBosonLep1CosThetaStarFJ = CosThetaStar(looseLeps->at(0)->p4(), looseLeps->at(1)->p4(), ZHP4);
    }
  }

  recluster->execute();
}


float FatJetMod::getMSDCorr(float puppipt, float puppieta) 
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

void FatJetReclusterMod::do_execute()
{
  if ((*fj1) == nullptr)
    return;

  VPseudoJet particles = convertPFCands(event.pfCandidates,analysis.puppiJets,0);
  fj::ClusterSequenceArea seq(particles,*jetDef,*(utils.areaDef));
  VPseudoJet allJets(seq.inclusive_jets(0.));
  fj::PseudoJet *pj1=0;
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


    fj::PseudoJet sdJet = (*utils.softDrop)(*pj1);
    VPseudoJet sdConstituents = fastjet::sorted_by_pt(sdJet.constituents());
    eTot=0; eTrunc=0;
    for (unsigned iC=0; iC!=sdConstituents.size(); ++iC) {
      double e = sdConstituents.at(iC).E();
      eTot += e;
      if (iC<100)
        eTrunc += e;
    }
  }
}

const GenParticle * FatJetMatchingMod::matchGen(double eta, double phi, double radius, int pdgid) const 
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


void FatJetMatchingMod::do_execute()
{
  auto* fj1 = *fjPtr;
  if (fj1 == nullptr)
    return;

  int pdgidTarget=0;
  if (!analysis.isData && analysis.processType>=kTT) {
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
        PError("FatJetMatchingMod::do_execute","Reached an unknown process type");
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

  if (!analysis.isData && gt.nFatjet>0) {
    // first see if jet is matched
    auto* matched = matchGen(fj1->eta(),fj1->phi(),1.5,pdgidTarget);
    if (matched!=nullptr) {
      gt.fjIsMatched = 1;
      gt.fjGenPt = matched->pt();
      gt.fjGenSize = genObjects[matched];
    } else {
      gt.fjIsMatched = 0;
    }
    if (pdgidTarget==6) { // matched to top; try for W
      auto* matchedW = matchGen(fj1->eta(),fj1->phi(),1.5,24);
      if (matchedW!=nullptr) {
        gt.fjIsWMatched = 1;
        gt.fjGenWPt = matchedW->pt();
        gt.fjGenWSize = genObjects[matchedW];
      } else {
        gt.fjIsWMatched = 0;
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
      if (pt>(gt.fjHighestPtGenPt)
          && DeltaR2(gen.eta(),gen.phi(),fj1->eta(),fj1->phi())<cfg.FATJETMATCHDR2) {
        gt.fjHighestPtGenPt = pt;
        gt.fjHighestPtGen = pdgid;
      }

      if (gen.parent.isValid() && gen.parent->pdgid==gen.pdgid)
        continue;

      //count bs and cs
      int apdgid = abs(pdgid);
      if (apdgid!=5 && apdgid!=4) 
        continue;

      if (DeltaR2(gen.eta(),gen.phi(),fj1->eta(),fj1->phi())<cfg.FATJETMATCHDR2) {
        gt.fjNHF++;
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

    gt.fjNbs=bs_inside_cone;
    gt.fjgbb=has_gluon_splitting;

    if (analysis.btagSFs) {
      // now get the subjet btag SFs
      vector<btagcand> sj_btagcands;
      vector<double> sj_sf_cent, sj_sf_bUp, sj_sf_bDown, sj_sf_mUp, sj_sf_mDown;
      int nSJ = fj1->subjets.size();
      for (int iSJ=0; iSJ!=nSJ; ++iSJ) {
        auto& subjet = fj1->subjets.objAt(iSJ);
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
