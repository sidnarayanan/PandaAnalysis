#include "../interface/BTagMods.h"

using namespace pa;
using namespace std;
using namespace panda;

void BTagSFMod::do_execute() 
{
  // now get the jet btag SFs
  vector<btagcand> btagcands;
  vector<double> sf_cent, sf_bUp, sf_bDown, sf_mUp, sf_mDown;

  auto& bcands = (*jesShifts)[0].bcand;
  auto& isojets = (*jesShifts)[0].iso;
  unsigned int nJ = bcands.size();

  for (unsigned int iJ=0; iJ!=nJ; ++iJ) {
    const JetWrapper *jw = bcands[iJ];
    const panda::Jet &jet = jw->get_base();

    bool isIsoJet = false;
    if (!analysis.fatjet || 
        // if we do not consider fatjets, everything is an isojet 
        // otherwise, explicitly check isojet
        find(isojets.begin(), isojets.end(), jw) != isojets.end()) 
      isIsoJet = true;

    int flavor = jw->flavor;
    float pt = jw->pt;
    float eta = jet.eta();

    float btagUncFactor = 1;
    double eff(1),sf(1),sfUp(1),sfDown(1);
    if (flavor==5)
      eff = utils.getCorr(cCSVBL,pt,fabs(eta));
    else if (flavor==4)
      eff = utils.getCorr(cCSVCL,pt,fabs(eta));
    else
      eff = utils.getCorr(cCSVLL,pt,fabs(eta));
    if (isIsoJet) {

      utils.btag->calcSF(bJetL,flavor,eta,pt,eff,btagUncFactor,sf,sfUp,sfDown);
      btagcands.emplace_back(iJ,flavor,eff,sf,sfUp,sfDown);
      sf_cent.push_back(sf);

      if (flavor>0) {
        sf_bUp.push_back(sfUp); sf_bDown.push_back(sfDown);
        sf_mUp.push_back(sf); sf_mDown.push_back(sf);
      } else {
        sf_bUp.push_back(sf); sf_bDown.push_back(sf);
        sf_mUp.push_back(sfUp); sf_mDown.push_back(sfDown);
      }

    }

  } // loop over jets

  utils.btag->evalSF(btagcands,sf_cent,GeneralTree::bCent,GeneralTree::bJet);
  utils.btag->evalSF(btagcands,sf_bUp,GeneralTree::bBUp,GeneralTree::bJet);
  utils.btag->evalSF(btagcands,sf_bDown,GeneralTree::bBDown,GeneralTree::bJet);
  utils.btag->evalSF(btagcands,sf_mUp,GeneralTree::bMUp,GeneralTree::bJet);
  utils.btag->evalSF(btagcands,sf_mDown,GeneralTree::bMDown,GeneralTree::bJet);
}

void BTagWeightMod::do_execute()
{
  auto& bcands = (*jesShifts)[0].bcand;
  if (bcands.size() < 1) 
    return;

  //get vectors of jet properties
  vector<double> jetPts, jetEtas, jetDiscrs;
  vector<int> jetFlavors;
  unsigned int nJ = bcands.size();
  jetPts.reserve(nJ);
  jetEtas.reserve(nJ);
  jetDiscrs.reserve(nJ);
  jetFlavors.reserve(nJ);
  for (unsigned int iJ=0; iJ!=nJ; ++iJ) {
    auto* jw = bcands[iJ];
    auto& jet = jw->get_base();
    jetPts.push_back(jw->pt);
    jetEtas.push_back(jet.eta());
    jetFlavors.push_back(jw->flavor);
    if (analysis.year==2016) {
      if (analysis.useCMVA)
        jetDiscrs.push_back(jet.cmva);
      else
        jetDiscrs.push_back(jet.csv);
    } else {
      jetDiscrs.push_back(jet.deepCSVb);
    }
  }
  for (unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
    GeneralTree::csvShift shift = gt.csvShifts[iShift];
    gt.sf_csvWeights[shift]=1;
    for (unsigned iJ=0; iJ!=nJ; ++iJ) {
      BTagEntry::JetFlavor flav;
      unsigned absid = abs(jetFlavors[iJ]);
      if     (absid==5) flav=BTagEntry::FLAV_B;
      else if(absid==4) flav=BTagEntry::FLAV_C;
      else              flav=BTagEntry::FLAV_UDSG;

      gt.sf_csvWeights[shift] *= utils.btag->reshaper->eval_auto_bounds(
        utils.btag->officialShiftNames[shift],
        flav, jetEtas[iJ], jetPts[iJ], jetDiscrs[iJ]);
    }
  }
}
