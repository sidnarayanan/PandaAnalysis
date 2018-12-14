#include "../interface/DeepOps.h"

using namespace pa;
using namespace std;
using namespace panda;
namespace fj = fastjet;
namespace tf = tensorflow;


float dnn_clean(float x) { 
  return fabs(x -  -99.) < 0.1 ? 0 : x; 
}

void BRegDeepOp::do_execute()
{
  auto& jw = **currentJet;
  int N = jw.cleaned_idx;

  // defined in SMH/breg/inputs.cfg
  inputs[ 0] = gt.jotRawPt[N];
  inputs[ 1] = gt.jotEta[N];
  inputs[ 2] = gt.jotPhi[N];
  inputs[ 3] = gt.jotRawMt[N];
  inputs[ 4] = gt.jotRawEt[N];
  inputs[ 5] = gt.jotRawM[N];
  inputs[ 6] = gt.jotRawE[N];
  inputs[ 7] = gt.jotPtD[N];
  inputs[ 8] = gt.jotCEF[N];
  inputs[ 9] = gt.jotNEF[N];
  inputs[10] = gt.jotCHF[N];
  inputs[11] = gt.jotNHF[N];
  inputs[12] = gt.jotNLep[N];
  inputs[13] = gt.jotLep1Pt[N];
  inputs[14] = gt.jotLep1Eta[N];
  inputs[15] = gt.jotLep1Phi[N];
  inputs[16] = gt.jotLep1PtRelRaw[N];
  inputs[17] = gt.jotLep1PtRelRawInv[N];
  inputs[18] = gt.jotLep1DeltaR[N];
  inputs[19] = gt.jotTrk1Pt[N];
  inputs[20] = gt.jotChTrk1Pt[N];
  inputs[21] = gt.jotVtxPt[N];
  inputs[22] = gt.jotVtxMass[N];
  inputs[23] = gt.jotVtx3DVal[N];
  inputs[24] = gt.jotVtx3DErr[N];
  inputs[25] = gt.jotVtxNtrk[N];
  inputs[26] = gt.jotLep1IsEle[N];
  inputs[27] = gt.jotLep1IsMu[N];
  inputs[28] = gt.jotLep1IsOther[N];
  inputs[29] = gt.jotNPt03[N];
  inputs[30] = gt.jotEMRing[0][N];
  inputs[31] = gt.jotEMRing[1][N];
  inputs[32] = gt.jotEMRing[2][N];
  inputs[33] = gt.jotEMRing[3][N];
  inputs[34] = gt.jotEMRing[4][N];
  inputs[35] = gt.jotChRing[0][N];
  inputs[36] = gt.jotChRing[1][N];
  inputs[37] = gt.jotChRing[2][N];
  inputs[38] = gt.jotChRing[3][N];
  inputs[39] = gt.jotChRing[4][N];
  inputs[40] = gt.jotMuRing[0][N];
  inputs[41] = gt.jotMuRing[1][N];
  inputs[42] = gt.jotMuRing[2][N];
  inputs[43] = gt.jotMuRing[3][N];
  inputs[44] = gt.jotMuRing[4][N];
  inputs[45] = gt.jotNeRing[0][N];
  inputs[46] = gt.jotNeRing[1][N];
  inputs[47] = gt.jotNeRing[2][N];
  inputs[48] = gt.jotNeRing[3][N];
  inputs[49] = gt.jotNeRing[4][N];
  inputs[50] = gt.rho;

  for (auto& i : inputs)
    i = dnn_clean(i); 

  eval();

  jw.breg = outputs[0];
  jw.bregwidth = 0.5 * (outputs[2] - outputs[1]);
}

void ZvvHClassOp::do_execute()
{
  // defined in SMH/evt/root.json
  inputs[ 0] = gt.hbbm_dreg[0];
  inputs[ 1] = abs(gt.jotEta[gt.hbbjtidx[0][1]]-gt.jotEta[gt.hbbjtidx[0][0]]);
  inputs[ 2] = gt.jotPt[0][gt.hbbjtidx[0][1]];
  inputs[ 3] = gt.jotPt[0][gt.hbbjtidx[0][0]];
  inputs[ 4] = gt.jotCSV[gt.hbbjtidx[0][1]];
  inputs[ 5] = gt.jotCSV[gt.hbbjtidx[0][0]];
  inputs[ 6] = fabs(SignedDeltaPhi(gt.jotPhi[gt.hbbjtidx[0][0]],
                                   gt.jotPhi[gt.hbbjtidx[0][1]]));
  inputs[ 7] = DeltaR2(gt.jotEta[gt.hbbjtidx[0][0]],
                       gt.jotPhi[gt.hbbjtidx[0][0]],
                       gt.jotEta[gt.hbbjtidx[0][1]],
                       gt.jotPhi[gt.hbbjtidx[0][1]]);
  inputs[ 8] = gt.nJot[0];
  inputs[ 9] = gt.nSoft5;
  inputs[10] = gt.adjetCMVA;
  inputs[11] = gt.adjetPt;

  eval();

  gt.zvvhClass = outputs[0];
}

void TFInferOp::build(TString weightpath)
{
  tf::setLogging("3");
  graph.reset(tf::loadGraphDef(weightpath.Data()));
  tf::SessionOptions opts; tf::setThreading(opts, 1, "no_threads");
  sess.reset(tf::createSession(graph.get(), opts));

  inputs.resize(n_inputs); 
  outputs.resize(n_outputs); 
  t_i[0] = tf::NamedTensor(inputName.Data(),
                           tf::Tensor(tf::DT_FLOAT,
                                      tf::TensorShape({1, (long long int)n_inputs})));
}

void TFInferOp::eval()
{
  // build the input tensor
  for (int idx = 0; idx != n_inputs; ++idx)
    t_i[0].second.matrix<float>()(0, idx) = inputs[idx];

  // set up output and run
  vector<tf::Tensor> t_o;
  tf::run(sess.get(), t_i, outputNames, &t_o);
  for (int idx = 0; idx != n_outputs; ++idx) {
    if (transpose)
      outputs[idx] = t_o[0].matrix<float>()(idx,0);
    else
      outputs[idx] = t_o[0].matrix<float>()(0,idx);
    //outputs[idx] = t_o[idx].tensor<float,2>()(0,0);
  }
}

template <typename GENP>
DeepGenOp<GENP>::DeepGenOp(panda::EventAnalysis& event_,
           Config& cfg_,
           Utils& utils_,
           GeneralTree& gt_,
           int level_) :
  AnalysisOp("deepgen", event_, cfg_, utils_, gt_, level_)
{
  if (!on())
    return;

  double radius = 1.5;
  auto algo = fj::cambridge_algorithm;
  if (analysis.ak8) {
    radius = 0.8;
    algo = fj::antikt_algorithm;
  }
  jetDef.reset(new fj::JetDefinition(algo,radius));
  tauN.reset(new fj::contrib::Njettiness(fj::contrib::OnePass_KT_Axes(),
                                         fj::contrib::NormalizedMeasure(1., radius)));
  ecfcalc.reset(new ECFCalculator());
  if (analysis.deepGenGrid) {
    grid.reset(new ParticleGridder(2500,1570,5)); // 0.002x0.002
    grid->_etaphi = false;
  }
}

template <typename GENP>
void DeepGenOp<GENP>::do_execute()
{
  gt.genFatJetPt = 0;
  gt.genFatJetNProngs = -1;

  set<int> leptonIndices; // hard leptons from t->W->lv ecay
  vector<fj::PseudoJet> finalStates;

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
  fj::ClusterSequenceArea seq(finalStates, *jetDef, *(utils.areaDef));
  vector<fj::PseudoJet> allJets(sorted_by_pt(seq.inclusive_jets(0.)));

  fj::PseudoJet* fullJet = NULL;
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
  fj::PseudoJet sdJet = (*utils.softDrop)(*fullJet);
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
  gt.fjTau32[0] = genJetInfo.tau3/genJetInfo.tau2;
  gt.fjTau21[0] = genJetInfo.tau2/genJetInfo.tau1;
  gt.fjTau32SD[0] = genJetInfo.tau3sd/genJetInfo.tau2sd;
  gt.fjTau21SD [0]= genJetInfo.tau2sd/genJetInfo.tau1sd;


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
    gt.fjECFNs[ep][0] = ecf;
  }

  // now we have to count the number of prongs
  unordered_set<const GENP*> partons;
  countGenPartons(partons); // fill the parton set

  map<const GENP*, unsigned> partonToIdx;
  for (auto* parton : partons)
    partonToIdx[parton] = partonToIdx.size(); // just some arbitrary ordering

  // get the hardest particle with angle wrt jet axis > 0.1
  fj::PseudoJet* axis2 = NULL;
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
    fj::PseudoJet &c = allConstituents.at(iC_);

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
void DeepGenOp<GENP>::countGenPartons(unordered_set<const GENP*>& partons)
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
void DeepGenOp<GENP>::incrementAux(bool close)
{
  if (fAux) {
    fAux->WriteTObject(tAux.get(), "inputs", "Overwrite");
    TString path = TString::Format(cfg.auxFilePath.Data(),auxCounter);
    fAux->Close();
  }
  if (close)
    return;

  TString path = TString::Format(cfg.auxFilePath.Data(),auxCounter++);
  fAux.reset(TFile::Open(path.Data(), "RECREATE"));
  tAux.reset(new TTree("inputs","inputs"));

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

template class DeepGenOp<panda::GenParticle>;
template class DeepGenOp<panda::UnpackedGenParticle>;

void BRegBDTOp::do_readData(TString dirPath)
{
  bjetreg_vars.reset(new float[15]);

  bjetregReader.reset(new TMVA::Reader("!Color:!Silent"));
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

  TString modelPath = dirPath + "/trainings/bjetregression.weights.xml";
  downloadData(
    "http://t3serv001.mit.edu/~dhsu/pandadata/trainings/bjet_regression_v1_fromBenedikt.weights.xml",
    modelPath.Data()
  );
  bjetregReader->BookMVA("BDT method", modelPath);
}

void BRegBDTOp::do_execute()
{
  auto& jw = **currentJet;
  int idx = jw.cleaned_idx;

  TLorentzVector v = jw.p4();

  // Shifted values for the jet energies to perform the b-jet regression
  bjetreg_vars[0] = jw.pt; // shifted pT
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
  bjetreg_vars[13]= v.Et();
  bjetreg_vars[14]= v.Mt();

  jw.breg = (bjetregReader->EvaluateRegression("BDT method"))[0];
  jw.bregwidth = 0;
}
