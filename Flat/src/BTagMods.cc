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

  for (unsigned iShift=0; iShift<GeneralTree::nCsvShifts; iShift++) {
    GeneralTree::csvShift shift = gt.csvShifts[iShift];
    gt.sf_csvWeights[shift]=1;
    for (auto* jw : bcands) {
      auto& jet = jw->get_base();
      float discr;
      if      (analysis.useCMVA   ) discr = jet.cmva;
      else if (analysis.useDeepCSV) discr = jet.deepCSVb + jet.deepCSVbb;
      else                          discr = jet.csv;
      unsigned absid = abs(jw->flavor);
      auto flav = absid == 5 ? BTagEntry::FLAV_B : 
                               (absid == 4 ? BTagEntry::FLAV_C : BTagEntry::FLAV_UDSG);
      float reshapeFactor = utils.btag->reshaper->eval_auto_bounds(
            GeneralTree::csvShiftName(shift).Data(),
            flav, jet.eta(), jw->pt, discr
            );
      if (reshapeFactor > 0.001)
        gt.sf_csvWeights[shift] *= reshapeFactor;
      else
        gt.sf_csvWeights[shift] *= utils.btag->reshaper->eval_auto_bounds(
            GeneralTree::csvShiftName(GeneralTree::csvCent).Data(),
            flav, jet.eta(), jw->pt, discr
            );
    }
  }
}
