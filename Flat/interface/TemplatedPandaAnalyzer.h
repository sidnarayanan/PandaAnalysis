#ifndef TEMPLATEPANDA
#define TEMPLATEPANDA 1 
    
template <typename T> 
void PandaAnalyzer::MatchGenJets(T& genJets) 
{
  JESLOOP {
    auto& jets = jesShifts[shift];
    unsigned N = jets.cleaned.size();
    for (unsigned i = 0; i != N; ++i) {
      const panda::Jet& reco = jets.cleaned[i]->get_base();
      for (auto &gen : genJets) {
        if (DeltaR2(gen.eta(), gen.phi(), reco.eta(), reco.phi()) < 0.09) {
          gt->jotGenPt[shift][i] = gen.pt();
          gt->jotFlav[shift][i] = gen.pdgid;
          break;
        }
      }
    }
  }
  tr->TriggerEvent("match gen jets");
}

template <typename T>
void PandaAnalyzer::RemoveGenDups(const panda::Collection<T>& genParticles)
{
  for (auto& g : genParticles) {
    bool foundDup = false;
    if (g.finalState) {
      float ptThreshold = g.pt() * 0.01;
      for (auto* pPtr : validGenP) {
        const T* gPtr = static_cast<const T*>(pPtr);
        if (!gPtr->finalState)
          continue;
        if ((g.pdgid == gPtr->pdgid) &&
            (fabs(g.pt() - gPtr->pt()) < ptThreshold) &&
            (DeltaR2(g.eta(), g.phi(), gPtr->eta(), gPtr->phi()) < 0.00001)) {
          foundDup = true;
          if (DEBUG > 8) {
            PDebug("Found duplicate",
                   Form("p1(%8.3f,%5.1f,%5.1f,%5i,%i) <-> p2(%8.3f,%5.1f,%5.1f,%5i,%i)",
                        g.pt(), g.eta(), g.phi(), g.pdgid, g.finalState ? 1 : 0,
                        gPtr->pt(), gPtr->eta(), gPtr->phi(), gPtr->pdgid, gPtr->finalState ? 1 : 0));
          }
          break;
        }
      }
    }
    if (!foundDup) {
      validGenP.push_back(&g);
    }
  } 
}

template <typename T>
void PandaAnalyzer::CountGenPartons(std::unordered_set<const T*>& partons)
{
  float dR2 = FATJETMATCHDR2;
  float base_eta = genJetInfo.eta, base_phi = genJetInfo.phi;
  auto matchJet = [base_eta, base_phi, dR2](const T& p) -> bool {
    return DeltaR2(base_eta, base_phi, p.eta(), p.phi()) < dR2;
  };
  float threshold = genJetInfo.pt * 0.2;
  for (auto* gen_ : validGenP) {
    auto* gen = dynamic_cast<const T*>(gen_);
    unsigned apdgid = abs(gen->pdgid);
    if (apdgid > 5 && 
        apdgid != 21 &&
        apdgid != 15 &&
        apdgid != 11 && 
        apdgid != 13)
      continue; 

    if (gen->pt() < threshold)
      continue; 

    if (!matchJet(*gen))
      continue;

    const T *parent = gen;
    const T *foundParent = NULL;
    while (parent->parent.isValid()) {
      parent = parent->parent.get();
      if (partons.find(parent) != partons.end()) {
        foundParent = parent;
        break;
      }
    }


    const T *dau1 = NULL, *dau2 = NULL;
    for (const auto* child_ : validGenP) {
      auto* child = dynamic_cast<const T*>(child_); 
      if (!(child->parent.isValid() && 
            child->parent.get() == gen))
        continue; 
      
      unsigned child_apdgid = abs(child->pdgid);
      if (child_apdgid > 5 && 
          child_apdgid != 21 &&
          child_apdgid != 15 &&
          child_apdgid != 11 && 
          child_apdgid != 13)
        continue; 

      if (dau1)
        dau2 = child;
      else
        dau1 = child;

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
      partons.insert(gen);
    }
  }

  genJetInfo.nprongs = partons.size();
  gt->genFatJetNProngs = genJetInfo.nprongs;

  TLorentzVector vPartonSum;
  TLorentzVector vTmp;
  for (auto *p : partons) {
    vTmp.SetPtEtaPhiM(p->pt(), p->eta(), p->phi(), p->m());
    vPartonSum += vTmp;
  }
  genJetInfo.partonm = vPartonSum.M();
  genJetInfo.partonpt = vPartonSum.Pt();
}


template <typename T>
void PandaAnalyzer::FillGenTree()
{

  genJetInfo.reset();  
  gt->genFatJetPt = 0;
  gt->genFatJetNProngs = -1;
  if (analysis->deepGenGrid)
    grid->clear();

  std::set<int> leptonIndices; // hard leptons from t->W->lv decay
  std::vector<fastjet::PseudoJet> finalStates;
  unsigned idx = -1;

  for (auto* p_ : validGenP) {
    auto* p = dynamic_cast<const T*>(p_);
    idx++;
    unsigned apdgid = abs(p->pdgid);
    if (!p->finalState)
      continue;
    if (apdgid == 12 ||
        apdgid == 14 ||
        apdgid == 16)
      continue; 
    if (p->pt() > 0.001 && fabs(p->eta()) < 5) {
      if (!analysis->deepGenGrid || (pdgToQ[apdgid] != 0)) { // it's charged, so we have tracking
        finalStates.emplace_back(p->px(), p->py(), p->pz(), p->e());
        finalStates.back().set_user_index(idx);

        if (apdgid == 11 ||
            apdgid == 13 ||
            apdgid == 15) {
          const T *parent = p;
          bool foundW = false, foundT = false;
          while (parent->parent.isValid()) {
            parent = parent->parent.get();
            unsigned parent_apdgid = abs(parent->pdgid);
            if (!foundW) {
              if (parent_apdgid == 24) {
                foundW = true;
                continue; 
              } else if (parent_apdgid != apdgid) {
                break; // if it's not a W, must be a parent of the particle we care about
              }
            } else {  // foundW = true
              if (parent_apdgid == 6) {
                foundT = true;
                break;
              } else if (parent_apdgid != 24) {
                break; // if it's not a top, must be a parent of the W we found
              }
            }
          }
          if (foundT) {
            leptonIndices.insert(idx);
          }
        }
      } else {
        grid->add(*p);
      }
    }
  }

  tr->TriggerSubEvent("reading gen");

  if (analysis->deepGenGrid) {
    int user_idx = -2;
    for (auto &v : grid->get()) {
      finalStates.emplace_back(v.Px(), v.Py(), v.Pz(), v.E());
      finalStates.back().set_user_index(user_idx); // not associated with a real particle
      --user_idx;
    } 
    tr->TriggerSubEvent("gen grid");
  }

  // to test exclusive clustering, don't need to cluster fat jets
  if (analysis->deepExC) {
    std::sort(finalStates.begin(), finalStates.end(),
              [](const fastjet::PseudoJet& x, const fastjet::PseudoJet& y) { return x.perp() > y.perp(); });
    for (int iC = 0; iC != std::min(NMAXPF,(int)finalStates.size()); ++iC) {
      auto& p = finalStates.at(iC);
      p.set_user_index(iC); // overwrite the indices
      genJetInfo.particles[iC][0] = p.eta();
      genJetInfo.particles[iC][1] = p.phi();
      genJetInfo.particles[iC][2] = p.perp();
      genJetInfo.particles[iC][3] = p.e();
    }
    // cluster into ak4 jets
    fastjet::ClusterSequenceArea seq(finalStates, *jetDefGen, *areaDef);
    std::vector<fastjet::PseudoJet> allJets(fastjet::sorted_by_pt(seq.inclusive_jets(0.)));
    for (int iJ = 0; iJ != std::min(5,(int)allJets.size()); ++iJ) {
      auto& j = allJets.at(iJ);
      if (j.perp() < 10)
        continue;
      for (auto& p : j.constituents()) {
        if (p.user_index() < NMAXPF && p.user_index() > 0) {
          genJetInfo.particles[p.user_index()][4] = iJ + 1;
        }
      }
    }
    tr->TriggerEvent("fill gen tree");
    return;
  }


  // cluster the  jet 
  fastjet::ClusterSequenceArea seq(finalStates, *jetDef, *areaDef);
  std::vector<fastjet::PseudoJet> allJets(fastjet::sorted_by_pt(seq.inclusive_jets(0.)));


  fastjet::PseudoJet* fullJet = NULL;
  for (auto& testJet : allJets) {
    if (testJet.perp() < minGenFatJetPt)
      break;
    bool jetOverlapsLepton = false;
    for (auto& c : testJet.constituents()) {
      int idx = c.user_index();
      if (leptonIndices.find(idx) != leptonIndices.end()) {
        jetOverlapsLepton = true;
        break;
      }
    }
    if (!jetOverlapsLepton) {
      fullJet = &testJet;
      break;
    }
  }

  tr->TriggerSubEvent("clustering gen");

  if (fullJet == NULL) {
    tr->TriggerEvent("fill gen tree");
    return;
  }

  gt->genFatJetPt = fullJet->perp();
  if (gt->genFatJetPt < 450) {
    tr->TriggerEvent("fill gen tree");
    return;
  }

  VPseudoJet allConstituents = fastjet::sorted_by_pt(fullJet->constituents());
  genJetInfo.pt = gt->genFatJetPt;
  genJetInfo.m = fullJet->m();
  genJetInfo.eta = fullJet->eta();
  genJetInfo.phi = fullJet->phi();

  // softdrop the jet
  fastjet::PseudoJet sdJet = (*softDrop)(*fullJet);
  VPseudoJet sdConstituents = fastjet::sorted_by_pt(sdJet.constituents());
  genJetInfo.msd = sdJet.m();
  std::vector<bool> survived(allConstituents.size());
  unsigned nC = allConstituents.size();
  for (unsigned iC = 0; iC != nC; ++iC) {
    int idx = allConstituents.at(iC).user_index();
    survived[iC] = false;
    for (auto& sdc : sdConstituents) {
      if (idx == sdc.user_index()) {
        survived[iC] = true; 
        break;
      }
    }
  }

  // get tau  
  genJetInfo.tau1 = tauN->getTau(1, allConstituents);
  genJetInfo.tau2 = tauN->getTau(2, allConstituents);
  genJetInfo.tau3 = tauN->getTau(3, allConstituents);
  genJetInfo.tau1sd = tauN->getTau(1, sdConstituents);
  genJetInfo.tau2sd = tauN->getTau(2, sdConstituents);
  genJetInfo.tau3sd = tauN->getTau(3, sdConstituents);
  gt->fjTau32 = genJetInfo.tau3/genJetInfo.tau2;
  gt->fjTau21 = genJetInfo.tau2/genJetInfo.tau1;
  gt->fjTau32SD = genJetInfo.tau3sd/genJetInfo.tau2sd;
  gt->fjTau21SD = genJetInfo.tau2sd/genJetInfo.tau1sd;

  tr->TriggerSubEvent("sd and tau");

  // get ecfs
  unsigned nFilter = std::min(30, (int)sdConstituents.size());
  VPseudoJet sdConstsFiltered(sdConstituents.begin(), sdConstituents.begin() + nFilter);

  GeneralTree::ECFParams ep;
  ecfcalc->calculate(sdConstsFiltered);
  for (auto iter = ecfcalc->begin(); iter != ecfcalc->end(); ++iter) {
    int N = iter.get<pandaecf::Calculator::nP>();
    int o = iter.get<pandaecf::Calculator::oP>();
    int beta = iter.get<pandaecf::Calculator::bP>();
    // float ecf = iter.get<pandaecf::Calculator::ecfP>().template convert_to<float>();
    float ecf(iter.get<pandaecf::Calculator::ecfP>());
    // PDebug("",Form("io=%i, iN=%i, ibeta=%i, ecf=%g", o, N, beta, ecf));
    genJetInfo.ecfs[o][N][beta] = ecf;
    ep.order = o + 1; ep.N = N + 1, ep.ibeta = beta;
    gt->fjECFNs[ep] = ecf;
  }
  
  tr->TriggerSubEvent("ecfs");


  // now we have to count the number of prongs 
  std::unordered_set<const T*> partons; 
  CountGenPartons<T>(partons); // fill the parton set

  std::map<const T*, unsigned> partonToIdx;
  for (auto* parton : partons) 
    partonToIdx[parton] = partonToIdx.size(); // just some arbitrary ordering 

  // get the hardest particle with angle wrt jet axis > 0.1
  fastjet::PseudoJet* axis2 = NULL;
  for (auto& c : allConstituents) {
    if (DeltaR2(c.eta(), c.phi(), genJetInfo.eta, genJetInfo.phi) > 0.01) {
      if (!axis2 || (c.perp() > axis2->perp())) {
        axis2 = &c;
      }
    }
  }

  JetRotation rot(fullJet->px(), fullJet->py(), fullJet->pz(),
                  axis2->px(), axis2->py(), axis2->pz());

  std::vector<int> indices;
  std::vector<int> unclustered; // I honestly don't know where these come from...
  if (analysis->deepAntiKtSort) {
    indices.reserve(nC);
    std::vector<bool> mask(nC, false);

    JetTree jt(*fullJet);
    tr->TriggerSubEvent("filling jet history");

    std::vector<int> terminalIdx = jt.GetTerminals();
    for (auto uidx : terminalIdx) {
      if (uidx != -1) {
        auto iter = std::find_if(allConstituents.begin(), allConstituents.end(),
                                 JetTree::compare(uidx));
        if (iter == allConstituents.end())
          continue;
        
        int idx = static_cast<int>(iter - allConstituents.begin());
        indices.push_back(idx);
        mask[idx] = true;
      }
    }

    unclustered.reserve(nC/5); // whatever, just an estimate
    for (unsigned iC = 0; iC != nC; ++iC) {
      if (!mask[iC])
        unclustered.push_back(iC);
    }

  }

  // now we fill the particles
  nC = std::min(nC, (unsigned)NMAXPF);
  for (unsigned iC = 0; iC != nC; ++iC) {
    unsigned iC_ = iC;
    if (analysis->deepAntiKtSort) {
      if (iC < indices.size()) {
        iC_ = indices.at(iC);
      } else {
        iC_ = unclustered.at(iC - indices.size());
      }
    }
    fastjet::PseudoJet &c = allConstituents.at(iC_);

    if (c.perp() < 0.001) // not a real particle
      continue;

    float angle = DeltaR2(c.eta(), c.phi(), genJetInfo.eta, genJetInfo.phi);
    float x=c.px(), y=c.py(), z=c.pz();
    rot.Rotate(x, y, z);  // perform two rotations on the jet 
    genJetInfo.particles[iC][0] = x;
    genJetInfo.particles[iC][1] = y;
    genJetInfo.particles[iC][2] = z;
    genJetInfo.particles[iC][3] = c.e();
    genJetInfo.particles[iC][4] = angle;
    genJetInfo.particles[iC][5] = survived[iC] ? 1 : 0;

    unsigned ptype = 0;
    int parent_idx = -1;
    if (c.user_index() >= 0) {
      const T* gen = dynamic_cast<const T*>(validGenP.at(c.user_index()));
      int pdgid = gen->pdgid;
      int apdgid = abs(pdgid);
      if (apdgid == 11) {
        ptype = 1 * sign(pdgid * -11);
      } else if (apdgid == 13) {
        ptype = 2 * sign(pdgid * -13);
      } else if (apdgid == 22) {
        ptype = 3;
      } else {
        float q = pdgToQ[apdgid];
        if (apdgid != pdgid)
          q *= -1;
        if (q == 0) 
          ptype = 4;
        else if (q > 0) 
          ptype = 5;
        else 
          ptype = 6;
      }

      const T *parent = gen;
      while (parent->parent.isValid()) {
        parent = parent->parent.get();
        if (partons.find(parent) != partons.end()) {
          parent_idx = partonToIdx[parent];
          break;
        }
      }
    }

    genJetInfo.particles[iC][6] = ptype;
    genJetInfo.particles[iC][7] = parent_idx;
  }

  tr->TriggerEvent("fill gen tree");
}


#endif
