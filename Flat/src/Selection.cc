#include "../interface/Selection.h"

static float best_recoil(const GeneralTree *gt, bool include_var) {
  float max_puppi = std::max({gt->puppimet, gt->puppiUZmag, gt->puppiUWmag, gt->puppiUAmag});
  float max_pf = std::max({gt->pfmet, gt->pfUZmag, gt->pfUWmag, gt->pfUAmag});
  if (!include_var)
    return std::max({max_puppi, max_pf});
  float max_pfUp = std::max({gt->pfmetUp, gt->pfUZmagUp, gt->pfUWmagUp, gt->pfUAmagUp});
  float max_pfDown = std::max({gt->pfmetDown, gt->pfUZmagDown, gt->pfUWmagDown, gt->pfUAmagDown});
  return std::max({max_puppi, max_pf, max_pfUp, max_pfDown});
}


bool LeptonSel::do_accept() const 
{ 
  if (gt->nLooseLep < 2)
    return false;
  if (gt->nLooseElectron >= 2 && gt->electronPt[0] > 20 && gt->electronPt[1] > 20)
    return true;
  if (gt->nLooseMuon >= 2 && gt->muonPt[0] > 20 && gt->muonPt[1] > 20)
    return true;
  if (gt->muonPt[0] > 20 && gt->muonPt[1] > 20)
    return true;
  return false;
}


bool LeptonFakeSel::do_accept() const 
{ 
  bool passFakeTrigger = (gt->trigger & (1<<kMuFakeTrig)) != 0 || (gt->trigger & (1<<kEleFakeTrig)) != 0;
  if (!passFakeTrigger)
    return false;
  return (gt->nLooseLep==1 ||
          (gt->nLooseLep==2 && gt->diLepMass > 70));
}


bool RecoilSel::do_accept() const
{
  float bu = best_recoil(gt, vary_jes);
  return (bu > threshold); 
}


bool MonotopSel::do_accept() const
{
  bool base = RecoilSel::do_accept();
  return (base && gt->nFatjet>0 && gt->fj1Pt>200);
}


bool MonohiggsSel::do_accept() const
{
  bool base = RecoilSel::do_accept();

  return (base && 
          (gt->hbbpt < 150 ||
           (gt->nFatjet>=1 && gt->fj1Pt>200)));
}


bool VHbbSel::do_accept() const
{
  float bestMet = std::max({gt->pfmetUp, gt->pfmetDown, gt->pfmet});
  float bestJet1 = std::max({gt->jet1PtUp, gt->jet1PtDown, gt->jet1Pt});
  float bestJet2 = std::max({gt->jet2PtUp, gt->jet2PtDown, gt->jet2Pt});

  // ZnnHbb
  if (bestMet>150 && bestJet1 > 50 && bestJet2 > 25 &&
      (gt->hbbpt>50 || gt->nFatjet>0)) 
  {
    return true;
  }

  // WlnHbb
  if (bestJet1>25 && bestJet2>25 &&
      ((gt->nTightElectron>0 && gt->electronPt[0]>25) ||
       (gt->nTightMuon>0 && gt->muonPt[0]>25)) &&
      (gt->hbbpt>50 || gt->nFatjet>0)) 
  {
    return true;
  }

  // ZllHbb
  if (bestJet1>25 && bestJet2>25 && 
      ((gt->nTightElectron>0 && gt->nLooseElectron>1 && gt->electronPt[0]>25 && gt->electronPt[1]>20) ||
       (gt->nTightMuon>0 && gt->nLooseMuon>1 && gt->muonPt[0]>25 && gt->muonPt[1]>20)) &&
      (gt->hbbpt>50 || gt->nFatjet>0))
  {
    return true;
  }

  return false;
}
