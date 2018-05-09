#include "../interface/DeepMods.h"

using namespace pa;
using namespace std;
using namespace panda;
using namespace fastjet;
namespace tf = tensorflow; 

void BRegDeepMod::do_execute()
{
  auto& jw = **currentJet; 
  auto& jet = jw.get_base();
  int N = jw.user_idx; 

  TLorentzVector vRaw;
  vRaw.SetPtEtaPhiM(jet.rawPt, jet.eta(), jet.phi(), jet.m());

  // defined in data/trainings/breg_training_2017.cfg
  inputs[ 0] = jet.rawPt; 
  inputs[ 1] = jet.eta();
  inputs[ 2] = event.rho;
  inputs[ 3] = vRaw.Mt();
  inputs[ 4] = gt.jotTrk1Pt[N];
  inputs[ 5] = gt.jotLep1PtRel[N];
  inputs[ 6] = gt.jotLep1DeltaR[N];
  inputs[ 7] = 0;
  inputs[ 8] = 0;
  inputs[ 9] = gt.jotVtxPt[N];
  inputs[10] = gt.jotVtxMass[N];
  inputs[11] = gt.jotVtx3DVal[N];
  inputs[12] = gt.jotVtxNtrk[N]; 
  inputs[13] = gt.jotVtx3DErr[N];
  inputs[14] = 0; 
  inputs[15] = 0;
  inputs[16] = 0;
  inputs[17] = 0;
  inputs[18] = 0;
  inputs[19] = 0;
  inputs[20] = 0;
  inputs[21] = 0;
  inputs[22] = 0;
  inputs[23] = 0;
  inputs[24] = 0;
  inputs[25] = 0;
  inputs[26] = 0;
  inputs[27] = 0;
  inputs[28] = 0;
  inputs[29] = 0;
  inputs[30] = 0;
  inputs[31] = 0;
  inputs[32] = 0;
  inputs[33] = 0;
  inputs[34] = 0;
  inputs[35] = 0;
  inputs[36] = 0;
  inputs[37] = 0;
  inputs[38] = 0;
  inputs[39] = 0;
  inputs[40] = 0;
  inputs[41] = 0;
  inputs[42] = jet.m(); 
  inputs[43] = 0;

  eval();

  jet.breg = outputs[0]*0.39077115058898926+1.0610932111740112;
  jet.bregwidth = 0.5*(outputs[2]-outputs[1])*0.39077115058898926;

}

void TFInferMod::build(TString weightpath)
{
  tf::setLogging("3");
  graph = tf::loadGraphDef(weightpath.Data());
  tf::SessionOptions opts; tf::setThreading(opts, 1, "no_threads");
  sess = tf::createSession(graph, opts); 

  inputs.resize(n_inputs); outputs.resize(n_outputs); 
  outputNames[0] = outputName.Data();
}

void TFInferMod::eval()
{
  // build the input tensor 
  // might be better to put this on the heap but need to check how it moves between sessions
  tf::NamedTensorList t_i; 
  t_i.resize(1);
  t_i[0] = tf::NamedTensor(inputName.Data(),
                           tf::Tensor(tf::DT_FLOAT, 
                                      tf::TensorShape({1, (long long int)n_inputs})));
  for (int idx = 0; idx != n_inputs; ++idx) 
    t_i[0].second.matrix<float>()(0, idx) = inputs[idx];

  // set up output and run
  vector<tf::Tensor> t_o; 
  tf::run(sess, t_i, outputNames, &t_o);
  for (int idx = 0; idx != n_outputs; ++idx)
    outputs[idx] = t_o[0].matrix<float>()(0, idx);
}

template <typename GENP>
DeepGenMod<GENP>::DeepGenMod(panda::EventAnalysis& event_, 
           Config& cfg_,
           Utils& utils_,
           GeneralTree& gt_) :
  AnalysisMod("deepgen", event_, cfg_, utils_, gt_) 
{ 
  if (!on())
    return;

  double radius = 1.5;
  auto algo = cambridge_algorithm;
  if (analysis.ak8) {
    radius = 0.8;
    algo = antikt_algorithm;
  } 
  jetDef = new JetDefinition(algo,radius);
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
  for (auto* p_ : *genP) {
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
          const GENP *parent = p;
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
    if (testJet.perp() < cfg.minGenFatJetPt)
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


  if (fullJet == NULL) {
    return;
  }

  gt.genFatJetPt = fullJet->perp();
  if (gt.genFatJetPt < 450) {
    return;
  }

  VPseudoJet allConstituents = sorted_by_pt(fullJet->constituents());
  genJetInfo.pt = gt.genFatJetPt;
  genJetInfo.m = fullJet->m();
  genJetInfo.eta = fullJet->eta();
  genJetInfo.phi = fullJet->phi();

  // softdrop the jet
  PseudoJet sdJet = (*utils.softDrop)(*fullJet);
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
    int N = iter.get<pa::ECFCalculator::nP>();
    int o = iter.get<pa::ECFCalculator::oP>();
    int beta = iter.get<pa::ECFCalculator::bP>();
    float ecf(iter.get<pa::ECFCalculator::ecfP>());
    genJetInfo.ecfs[o][N][beta] = ecf;
    ep.order = o + 1; ep.N = N + 1, ep.ibeta = beta;
    gt.fjECFNs[ep] = ecf;
  }
  
  // now we have to count the number of prongs 
  unordered_set<const GENP*> partons; 
  countGenPartons(partons); // fill the parton set

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
      auto* gen = dynamic_cast<const GENP*>(genP->at(c.user_index()));
      int pdgid = gen->pdgid;
      int apdgid = abs(pdgid);
      if (apdgid == 11) {
        ptype = -1 * sign(pdgid);
      } else if (apdgid == 13) {
        ptype = -2 * sign(pdgid);
      } else if (apdgid == 22) {
        ptype = 3;
      } else {
        float q = utils.pdgToQ[apdgid];
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
  for (auto* gen_ : *genP) {
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

    const GENP *parent = gen;
    const GENP *foundParent = NULL;
    while (parent->parent.isValid()) {
      parent = parent->parent.get();
      if (partons.find(parent) != partons.end()) {
        foundParent = parent;
        break;
      }
    }

    const GENP *dau1 = NULL, *dau2 = NULL;
    for (const auto* child_ : *genP) {
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
    TString path = TString::Format(cfg.auxFilePath.Data(),auxCounter);
    fAux->Close();
  }
  if (close)
    return;

  TString path = TString::Format(cfg.auxFilePath.Data(),auxCounter++);
  fAux = TFile::Open(path.Data(), "RECREATE");
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

template class DeepGenMod<panda::GenParticle>;
template class DeepGenMod<panda::UnpackedGenParticle>;
