#include "../interface/PandaAnalyzer.h"
#include "TVector2.h"
#include "TMath.h"
#include <algorithm>
#include <vector>
#include <unordered_set>


using namespace panda;
using namespace std;

// determine the partons associated with a fatjet
// Responsible: S. Narayanan
void PandaAnalyzer::FatjetPartons() 
{
  gt->fjNPartons = 0;
  if (fj1) {
    double threshold = 0.2 * fj1->rawPt;

    const FatJet *my_fj = fj1; 
    double dR2 = FATJETMATCHDR2; // put these guys in local scope

    auto matchJet = [my_fj, dR2](const GenParticle &p) -> bool {
      return DeltaR2(my_fj->eta(), my_fj->phi(), p.eta(), p.phi()) < dR2;
    };

    unordered_set<const panda::GenParticle*> partons; 
    for (auto* genptr : validGenP) {
      auto& gen = pToGRef(genptr);
      int apdgid = abs(gen.pdgid);
      if (apdgid > 5 && 
          apdgid != 21 &&
          apdgid != 15 &&
          apdgid != 11 && 
          apdgid != 13)
        continue; 

      if (gen.pt() < threshold)
        continue; 

      if (!matchJet(gen))
        continue;

      const GenParticle *parent = &gen;
      const GenParticle *foundParent = nullptr;
      while (parent->parent.isValid()) {
        parent = parent->parent.get();
        if (partons.find(parent) != partons.end()) {
          foundParent = parent;
          break;
        }
      }


      const GenParticle *dau1 = nullptr, *dau2 = nullptr;
      for (auto* childptr : validGenP) {
        auto& child = pToGRef(childptr);
        if (!(child.parent.isValid() && 
              child.parent.get() == &gen))
          continue; 
        
        int child_apdgid = abs(child.pdgid);
        if (child_apdgid > 5 && 
            child_apdgid != 21 &&
            child_apdgid != 15 &&
            child_apdgid != 11 && 
            child_apdgid != 13)
          continue; 

        if (dau1)
          dau2 = &child;
        else
          dau1 = &child;

        if (dau1 && dau2)
          break;
      }

      if (dau1 && dau2 && 
          dau1->pt() > threshold && dau2->pt() > threshold && 
          matchJet(*dau1) && matchJet(*dau2)) {
        if (foundParent) {
          partons.erase(partons.find(foundParent));
        }
        partons.insert(dau1);
        partons.insert(dau2);
      } else if (foundParent) {
        continue; 
      } else {
        partons.insert(&gen);
      }
    }

    gt->fjNPartons = partons.size();

    TLorentzVector vPartonSum;
    TLorentzVector vTmp;
    for (auto *p : partons) {
      vTmp.SetPtEtaPhiM(p->pt(), p->eta(), p->phi(), p->m());
      vPartonSum += vTmp;

      int digit3 = (p->pdgid%1000 - p->pdgid%100) / 100;
      int digit4 = (p->pdgid%10000 - p->pdgid%1000) / 1000;
      if (p->pdgid == 5 || digit3 == 5 || digit4 == 5)
        gt->fjNBPartons++;
      if (p->pdgid == 4 || digit3 == 4 || digit4 == 4)
        gt->fjNCPartons++;
    }
    gt->fjPartonM = vPartonSum.M();
    gt->fjPartonPt = vPartonSum.Pt();
    gt->fjPartonEta = vPartonSum.Eta();
  }

}

// fill an aux tree with reco info of the fatjet
// Responsible: S. Narayanan
void PandaAnalyzer::FillPFTree() 
{
  // this function saves the PF information of the leading fatjet
  // to a compact 2D array, that is eventually saved to an auxillary
  // tree/file. this is used as temporary input to the inference step
  // which then adds a separate tree to the main output. ideally,
  // this is integrated by use of e.g. lwtnn, but there is a bit of
  // development needed to support our networks. for the time being
  // we do it this way. -SMN
   

  // jet-wide quantities
  fjpt = -1; fjmsd = -1; fjeta = -1; fjphi = -1; fjphi = -1; fjrawpt = -1;
  for (int i = 0; i != NMAXPF; ++i) {
    for (int j = 0; j != NPFPROPS; ++j) {
      pfInfo[i][j] = 0;
    }
  }
  if (analysis->deepSVs) {
    for (int i = 0; i != NMAXSV; ++i) {
      for (int j = 0; j != NSVPROPS; ++j) {
        svInfo[i][j] = 0;
      }
    }}

  if (!fj1)
    return;

  fjpt = fj1->pt();
  fjmsd = fj1->mSD;
  fjeta = fj1->eta();
  fjphi = fj1->phi();
  fjrawpt = fj1->rawPt;

  gt->fjRho = TMath::Log(TMath::Power(fjmsd,2) / fjpt);
  gt->fjRawRho = TMath::Log(TMath::Power(fjmsd,2) / fjrawpt);
  gt->fjRho2 = TMath::Log(TMath::Power(fjmsd,2) / TMath::Power(fjpt,2));
  gt->fjRawRho2 = TMath::Log(TMath::Power(fjmsd,2) / TMath::Power(fjrawpt,2));

  // secondary vertices come first so we can link to tracks later
  std::map<const SecondaryVertex*, std::unordered_set<const PFCand*>> svTracks;
  std::map<const SecondaryVertex*, int> svIdx;
  if (analysis->deepSVs) {
    int idx = 0;
    for (auto &sv : event.secondaryVertices) {
      if (idx == NMAXSV)
        break;

      svInfo[idx][0] = sv.x;
      svInfo[idx][1] = sv.y;
      svInfo[idx][2] = sv.z;
      svInfo[idx][3] = sv.ntrk;
      svInfo[idx][4] = sv.ndof;
      svInfo[idx][5] = sv.chi2;
      svInfo[idx][6] = sv.significance;
      svInfo[idx][7] = sv.vtx3DVal;
      svInfo[idx][8] = sv.vtx3DeVal;
      svInfo[idx][9] = sv.pt();
      svInfo[idx][10] = sv.eta();
      svInfo[idx][11] = sv.phi();
      svInfo[idx][12] = sv.m();

      auto *svPtr = &sv;
      svIdx[svPtr] = idx;
      svTracks[svPtr] = {};
      for (auto pf : sv.daughters) {
        svTracks[svPtr].insert(pf.get());
      }
    }
  }


  std::vector<const PFCand*> sortedC;
  if (analysis->deepKtSort || analysis->deepAntiKtSort) {
    VPseudoJet particles = ConvertPFCands(fj1->constituents,analysis->puppiJets,0.001);
    fastjet::ClusterSequenceArea seq(particles,
                                     *(analysis->deepKtSort ? jetDefKt : jetDef),
                                     *areaDef);
    VPseudoJet allJets(seq.inclusive_jets(0.));
  
    auto &history = seq.history();
    auto &jets = seq.jets();
    vector<JetHistory> ordered_jets;
    for (auto &h : history) {
      if (h.jetp_index >= 0) {
        auto &j = jets.at(h.jetp_index);
        if (j.user_index() >= 0) {
          JetHistory jh;
          jh.user_idx = j.user_index();
          jh.child_idx = h.child;
          ordered_jets.push_back(jh);
        }
      }
    }
    sort(ordered_jets.begin(), ordered_jets.end(),
         [](JetHistory x, JetHistory y) { return x.child_idx < y.child_idx; });
    for (auto &jh : ordered_jets) {
      const PFCand *cand = fj1->constituents.at(jh.user_idx).get();
      sortedC.push_back(cand);
    }
  } else {
    for (auto ref : fj1->constituents)
      sortedC.push_back(ref.get());
    sort(sortedC.begin(), sortedC.end(),
         [](const PFCand *x, const PFCand *y) { return x->pt() > y->pt(); });
  }

  int idx = 0;
  for (auto *cand : sortedC) {
    if (idx == NMAXPF)
      break;
    pfInfo[idx][0] = cand->pt() * cand->puppiW() / fjrawpt;
    pfInfo[idx][1] = cand->eta() - fj1->eta();
    pfInfo[idx][2] = SignedDeltaPhi(cand->phi(), fj1->phi());
    pfInfo[idx][3] = cand->m();
    pfInfo[idx][4] = cand->e();
    pfInfo[idx][5] = cand->ptype;
    pfInfo[idx][6] = cand->puppiW();
    pfInfo[idx][7] = cand->puppiWNoLep(); 
    pfInfo[idx][8] = cand->hCalFrac;
    if (analysis->deepTracks) {
      if (cand->track.isValid())  {
        TVector3 pca = cand->pca();
        pfInfo[idx][9] = pca.Perp();
        pfInfo[idx][10] = pca.Z();
        pfInfo[idx][11] = cand->track->ptError();
        pfInfo[idx][12] = cand->track->dxy();
        pfInfo[idx][13] = cand->track->dz();
        pfInfo[idx][14] = cand->track->dPhi();
        pfInfo[idx][15] = cand->q();
        if (analysis->deepSVs) {
          for (auto &iter : svTracks) {
            if (iter.second.find(cand) != iter.second.end()) {
              TVector3 pos = iter.first->position();
              pfInfo[idx][16] = svIdx[iter.first] + 1; // offset from 0 value
              pfInfo[idx][17] = cand->dxy(pos);
              pfInfo[idx][18] = cand->dz(pos); 
              break;
            }
          }
        }
      }
    }
    idx++;
  }

  tr->TriggerEvent("pf tree");

}



// identify interesting gen particles for fatjet matching
// Responsible: S. Narayanan
void PandaAnalyzer::FatjetMatching() 
{
  int pdgidTarget=0;
  if (!isData && analysis->processType>=kTT) {
    switch(analysis->processType) {
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
        // analysis->processType>=kTT means we should never get here
        PError("PandaAnalyzer::Run","Reached an unknown process type");
    }

    std::vector<int> targets;

    int nGen = validGenP.size();
    for (int iG=0; iG!=nGen; ++iG) {
      auto& part = pToGRef(validGenP.at(iG));
      int pdgid = part.pdgid;
      int abspdgid = abs(pdgid);
      if (abspdgid == pdgidTarget)
        targets.push_back(iG);
    } //looking for targets

    for (int iG : targets) {
      auto& part = pToGRef(validGenP.at(iG));

      // check there is no further copy:
      bool isLastCopy=true;
      for (int jG : targets) {
        if (pToGPtr(validGenP.at(jG))->parent.get() == &part) {
          isLastCopy=false;
          break;
        }
      }
      if (!isLastCopy)
        continue;

      // (a) check it is a hadronic decay and if so, (b) calculate the size
      if (analysis->processType==kTop||analysis->processType==kTT) {

        // first look for a W whose parent is the top at iG, or a W further down the chain
        panda::GenParticle const* lastW(0);
        for (int jG=0; jG!=nGen; ++jG) {
          GenParticle const& partW = pToGRef(validGenP.at(jG));
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
          auto& partQ = pToGRef(validGenP.at(jG));
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
          auto& partQ = pToGRef(validGenP.at(jG));
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

  tr->TriggerEvent("gen matching");

  if (!isData && gt->nFatjet>0) {
    // first see if jet is matched
    auto* matched = MatchToGen(fj1->eta(),fj1->phi(),1.5,pdgidTarget);
    if (matched!=nullptr) {
      gt->fjIsMatched = 1;
      gt->fjGenPt = matched->pt();
      gt->fjGenSize = genObjects[matched];
    } else {
      gt->fjIsMatched = 0;
    }
    if (pdgidTarget==6) { // matched to top; try for W
      auto* matchedW = MatchToGen(fj1->eta(),fj1->phi(),1.5,24);
      if (matchedW!=nullptr) {
        gt->fjIsWMatched = 1;
        gt->fjGenWPt = matchedW->pt();
        gt->fjGenWSize = genObjects[matchedW];
      } else {
        gt->fjIsWMatched = 0;
      }
    }

    bool found_b_from_g=false;
    int bs_inside_cone=0;
    int has_gluon_splitting=0;
    panda::GenParticle const* first_b_mo(0);
    // now get the highest pT gen particle inside the jet cone
    for (auto* genptr : validGenP) {
      auto& gen = pToGRef(genptr);
      float pt = gen.pt();
      int pdgid = gen.pdgid;
      if (pt>(gt->fjHighestPtGenPt)
          && DeltaR2(gen.eta(),gen.phi(),fj1->eta(),fj1->phi())<FATJETMATCHDR2) {
        gt->fjHighestPtGenPt = pt;
        gt->fjHighestPtGen = pdgid;
      }

      if (gen.parent.isValid() && gen.parent->pdgid==gen.pdgid)
        continue;

      //count bs and cs
      int apdgid = abs(pdgid);
      if (apdgid!=5 && apdgid!=4) 
        continue;

      if (DeltaR2(gen.eta(),gen.phi(),fj1->eta(),fj1->phi())<FATJETMATCHDR2) {
        gt->fjNHF++;
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

    gt->fjNbs=bs_inside_cone;
    gt->fjgbb=has_gluon_splitting;

    if (analysis->btagSFs) {
      // now get the subjet btag SFs
      vector<btagcand> sj_btagcands;
      vector<double> sj_sf_cent, sj_sf_bUp, sj_sf_bDown, sj_sf_mUp, sj_sf_mDown;
      int nSJ = fj1->subjets.size();
      for (int iSJ=0; iSJ!=nSJ; ++iSJ) {
        auto& subjet = fj1->subjets.objAt(iSJ);
        int flavor=0;
        for (auto* genptr : validGenP) {
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
        int binpt = btagpt.bin(pt);
        int bineta = btageta.bin(fabs(eta));
        if (flavor==5) {
          eff = beff[bineta][binpt];
        } else if (flavor==4) {
          eff = ceff[bineta][binpt];
        } else {
          eff = lfeff[bineta][binpt];
        }
        if (analysis->hbb)
          CalcBJetSFs(bSubJetM,flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
        else
          CalcBJetSFs(bSubJetL,flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
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
      EvalBTagSF(sj_btagcands,sj_sf_cent,GeneralTree::bCent,GeneralTree::bSubJet,true);
      EvalBTagSF(sj_btagcands,sj_sf_bUp,GeneralTree::bBUp,GeneralTree::bSubJet,true);
      EvalBTagSF(sj_btagcands,sj_sf_bDown,GeneralTree::bBDown,GeneralTree::bSubJet,true);
      EvalBTagSF(sj_btagcands,sj_sf_mUp,GeneralTree::bMUp,GeneralTree::bSubJet,true);
      EvalBTagSF(sj_btagcands,sj_sf_mDown,GeneralTree::bMDown,GeneralTree::bSubJet,true);
    }

  }

  tr->TriggerEvent("fatjet gen-matching");
}

