#include "../interface/DeepMods.h"

using namespace pa;
using namespace std;
using namespace panda;
using namespace fastjet;

typedef DeepGenMod<GenParticle> DeepPGenMod;
typedef DeepGenMod<UnpackedGenParticle> DeepUGenMod;

template <typename GENP>
DeepGenMod<GENP>::DeepGenMod(panda::EventAnalysis& event_, 
           const Config& cfg_,
           const Utils& utils_,
           GeneralTree& gt_) :
  AnalysisMod("deepgen", event_, cfg_, utils_, gt_) 
{ 
  if (!on())
    return;

  double radius = 1.5;
  double sdZcut = 0.15;
  double sdBeta = 1.;
  auto algo = cambridge_algorithm;
  if (analysis.ak8) {
    radius = 0.8;
    sdZcut = 0.1;
    sdBeta = 0.;
    algo = antikt_algorithm;
  } 
  jetDef = new JetDefinition(algo,radius);
  softDrop = new contrib::SoftDrop(sdBeta,sdZcut,radius);
  tauN = new contrib::Njettiness(contrib::OnePass_KT_Axes(), 
                                          contrib::NormalizedMeasure(1., radius));
  ecfcalc = new ECFCalculator();
  if (analysis.deepGenGrid) {
    grid = new ParticleGridder(2500,1570,5); // 0.002x0.002
    grid->_etaphi = false;
  }
}

template <typename GENP>
void DeepGenMod<GENP>::do_execute() 
{
  gt.genFatJetPt = 0;
  gt.genFatJetNProngs = -1;

  set<int> leptonIndices; // hard leptons from t->W->lv ecay
  vector<PseudoJet> finalStates;

  int idx = -1;
  for (auto* p_ : genP) {
    auto* p = dynamic_cast<const GENP*>(p_);
    ++idx;
    if (!p->finalState)
      continue;
    unsigned apdgid = abs(p->pdgid);
    if (apdgid == 12 ||
        apdgid == 14 ||
        apdgid == 16)
      continue; 
    if (p->pt() > 0.001 && fabs(p->eta()) < 5) {
      if (!analysis.deepGenGrid || (utils.pdgToQ[apdgid] != 0)) { // it's charged, so we have tracking
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

  if (analysis.deepGenGrid) {
    int user_idx = -2;
    for (auto &v : grid->get()) {
      finalStates.emplace_back(v.Px(), v.Py(), v.Pz(), v.E());
      finalStates.back().set_user_index(user_idx); // not associated with a real particle
      --user_idx;
    } 
  }

  // cluster the  jet 
  ClusterSequenceArea seq(finalStates, *jetDef, *(utils.areaDef));
  vector<PseudoJet> allJets(sorted_by_pt(seq.inclusive_jets(0.)));

  PseudoJet* fullJet = NULL;
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

  gt.genFatJetPt = fullJet->perp();
  if (gt.genFatJetPt < 450) {
    tr->TriggerEvent("fill gen tree");
    return;
  }

  VPseudoJet allConstituents = sorted_by_pt(fullJet->constituents());
  genJetInfo.pt = gt.genFatJetPt;
  genJetInfo.m = fullJet->m();
  genJetInfo.eta = fullJet->eta();
  genJetInfo.phi = fullJet->phi();

  // softdrop the jet
  PseudoJet sdJet = (*softDrop)(*fullJet);
  VPseudoJet sdConstituents = sorted_by_pt(sdJet.constituents());
  genJetInfo.msd = sdJet.m();
  vector<bool> survived(allConstituents.size());
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
  gt.fjTau32 = genJetInfo.tau3/genJetInfo.tau2;
  gt.fjTau21 = genJetInfo.tau2/genJetInfo.tau1;
  gt.fjTau32SD = genJetInfo.tau3sd/genJetInfo.tau2sd;
  gt.fjTau21SD = genJetInfo.tau2sd/genJetInfo.tau1sd;


  // get ecfs
  unsigned nFilter = min(30, (int)sdConstituents.size());
  VPseudoJet sdConstsFiltered(sdConstituents.begin(), sdConstituents.begin() + nFilter);

  GeneralTree::ECFParams ep;
  ecfcalc->calculate(sdConstsFiltered);
  for (auto iter = ecfcalc->begin(); iter != ecfcalc->end(); ++iter) {
    int N = iter.get<pandaecf::Calculator::nP>();
    int o = iter.get<pandaecf::Calculator::oP>();
    int beta = iter.get<pandaecf::Calculator::bP>();
    float ecf(iter.get<pandaecf::Calculator::ecfP>());
    genJetInfo.ecfs[o][N][beta] = ecf;
    ep.order = o + 1; ep.N = N + 1, ep.ibeta = beta;
    gt.fjECFNs[ep] = ecf;
  }
  
  // now we have to count the number of prongs 
  unordered_set<const GENP*> partons; 
  countGenPartons<GENP>(partons); // fill the parton set

  map<const GENP*, unsigned> partonToIdx;
  for (auto* parton : partons) 
    partonToIdx[parton] = partonToIdx.size(); // just some arbitrary ordering 

  // get the hardest particle with angle wrt jet axis > 0.1
  PseudoJet* axis2 = NULL;
  for (auto& c : allConstituents) {
    if (DeltaR2(c.eta(), c.phi(), genJetInfo.eta, genJetInfo.phi) > 0.01) {
      if (!axis2 || (c.perp() > axis2->perp())) {
        axis2 = &c;
      }
    }
  }

  JetRotation rot(fullJet->px(), fullJet->py(), fullJet->pz(),
                  axis2->px(), axis2->py(), axis2->pz());

  vector<int> indices;
  vector<int> unclustered; // I honestly don't know where these come from...
  if (analysis.deepAntiKtSort) {
    indices.reserve(nC);
    vector<bool> mask(nC, false);

    JetTree jt(*fullJet);

    vector<int> terminalIdx = jt.GetTerminals();
    for (auto uidx : terminalIdx) {
      if (uidx != -1) {
        auto iter = find_if(allConstituents.begin(), allConstituents.end(),
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
  nC = min(nC, (unsigned)cfg.NMAXPF);
  for (unsigned iC = 0; iC != nC; ++iC) {
    unsigned iC_ = iC;
    if (analysis.deepAntiKtSort) {
      if (iC < indices.size()) {
        iC_ = indices.at(iC);
      } else {
        iC_ = unclustered.at(iC - indices.size());
      }
    }
    PseudoJet &c = allConstituents.at(iC_);

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
      auto* gen = dynamic_cast<const GENP*>(validGenP.at(c.user_index()));
      int pdgid = gen->pdgid;
      int apdgid = abs(pdgid);
      if (apdgid == 11) {
        ptype = -1 * sign(pdgid);
      } else if (apdgid == 13) {
        ptype = -2 * sign(pdgid);
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

      const GENP* parent = gen;
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

  if (gt.genFatJetPt > 450) {
    tAux->Fill();
    if (tAux->GetEntries() == 2500) 
      incrementAux(false);
  }
}

template <typename GENP>
void DeepGenMod<GENP>::countGenPartons(unordered_set<const GENP*>& partons) 
{
  float dR2 = cfg.FATJETMATCHDR2;
  float base_eta = genJetInfo.eta, base_phi = genJetInfo.phi;
  auto matchJet = [base_eta, base_phi, dR2](const GENP& p) -> bool {
    return DeltaR2(base_eta, base_phi, p.eta(), p.phi()) < dR2;
  };
  float threshold = genJetInfo.pt * 0.2;
  for (auto* gen_ : genP) {
    auto* gen = dynamic_cast<const GENP*>(gen_);
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
      auto* child = dynamic_cast<const GENP*>(child_); 
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
  gt.genFatJetNProngs = genJetInfo.nprongs;

  TLorentzVector vPartonSum;
  TLorentzVector vTmp;
  for (auto *p : partons) {
    vTmp.SetPtEtaPhiM(p->pt(), p->eta(), p->phi(), p->m());
    vPartonSum += vTmp;
  }
  genJetInfo.partonm = vPartonSum.M();
  genJetInfo.partonpt = vPartonSum.Pt();
}

template <typename GENP>
void DeepGenMod<GENP>::incrementAux(bool close)
{
  if (fAux) {
    fAux->WriteTObject(tAux, "inputs", "Overwrite");
    TString path = TString::Format(auxFilePath.Data(),auxCounter);
    if (DEBUG) PDebug("DeepGenMod::incrementAux", "Closing "+path);
    fAux->Close();
  }
  if (close)
    return;

  TString path = TString::Format(auxFilePath.Data(),auxCounter++);
  fAux = TFile::Open(path.Data(), "RECREATE");
  if (DEBUG) PDebug("DeepGenMod::incrementAux", "Opening "+path);
  tAux = new TTree("inputs","inputs");
  
  genJetInfo.particles.resize(cfg.NMAXPF);
  for (int i = 0; i != cfg.NMAXPF; ++i) {
    genJetInfo.particles[i].resize(cfg.NGENPROPS);
  }

  genJetInfo.ecfs.resize(3);
  for (int o = 0; o != 3; ++o) {
    genJetInfo.ecfs[o].resize(4);
    for (int N = 0; N != 4; ++N) {
      genJetInfo.ecfs[o][N].resize(4);
    }
  }

  tAux->Branch("eventNumber",&(gt.eventNumber),"eventNumber/l");
  tAux->Branch("nprongs",&(genJetInfo.nprongs),"nprongs/I");
  tAux->Branch("partonpt",&(genJetInfo.partonpt),"partonpt/F");
  tAux->Branch("partonm",&(genJetInfo.partonm),"partonm/F");
  tAux->Branch("pt",&(genJetInfo.pt),"pt/F");
  tAux->Branch("msd",&(genJetInfo.msd),"msd/F");
  tAux->Branch("eta",&(genJetInfo.eta),"eta/F");
  tAux->Branch("phi",&(genJetInfo.phi),"phi/F");
  tAux->Branch("m",&(genJetInfo.m),"m/F");
  tAux->Branch("tau3",&(genJetInfo.tau3),"tau3/F");
  tAux->Branch("tau2",&(genJetInfo.tau2),"tau2/F");
  tAux->Branch("tau1",&(genJetInfo.tau1),"tau1/F");
  tAux->Branch("tau3sd",&(genJetInfo.tau3sd),"tau3sd/F");
  tAux->Branch("tau2sd",&(genJetInfo.tau2sd),"tau2sd/F");
  tAux->Branch("tau1sd",&(genJetInfo.tau1sd),"tau1sd/F");
  for (int o = 1; o != 4; ++o) {
    for (int N = 1; N != 5; ++N) {
      for (int beta = 1; beta != 3; ++beta) {
        TString bname = Form("%i_%i_%i",o,N,beta);
        tAux->Branch(bname,&(genJetInfo.ecfs[o-1][N-1][beta-1]),bname+"/F");
      }
    }
  }
  tAux->Branch("kinematics",&(genJetInfo.particles));

  fOut->cd();
}
