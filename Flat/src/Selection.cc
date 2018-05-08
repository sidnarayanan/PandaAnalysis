#include "../interface/Selection.h"
#include "../interface/Common.h"

using namespace pa; 

static float best_recoil(const GeneralTree *gt, bool include_var) {
  int maxshift = include_var ? static_cast<int>(shiftjes::N) : 1;
  float max_puppi = -1, max_pf = -1;
  for (int shift = 0; shift != maxshift; ++ shift) {
    max_puppi = std::max({max_puppi, gt->puppimet[shift], 
                          gt->puppiUZmag[shift], 
                          gt->puppiUWmag[shift], 
                          gt->puppiUAmag[shift]});
    max_pf = std::max({max_pf, gt->pfmet[shift], 
                          gt->pfUZmag[shift], 
                          gt->pfUWmag[shift], 
                          gt->pfUAmag[shift]});
  }
  return std::max(max_puppi, max_pf);
}


bool LeptonSel::do_accept() const 
{ 
  if (gt->nLooseLep < 2)
    return false;
  if (gt->nLooseElectron >= 2 && gt->electronPt[0] > 20 && gt->electronPt[1] > 20)
    return true;
  if (gt->nLooseMuon >= 2 && gt->muonPt[0] > 20 && gt->muonPt[1] > 20)
    return true;
  if (gt->nLooseMuon >= 1 && gt->nLooseElectron >= 1 && gt->muonPt[0] > 20 && gt->electronPt[0] > 20)
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
  return (base && gt->nFatjet>0 && gt->fjPt[0]>200);
}


bool MonohiggsSel::do_accept() const
{
  bool base = RecoilSel::do_accept();

  return (base && 
          (gt->hbbpt[0] < 150 ||
           (gt->nFatjet>=1 && gt->fjPt[0]>200)));
}


bool VHbbSel::do_accept() const
{
  int maxshift = static_cast<int>(shiftjes::N);
  float bestMet = -1, bestJet1 = -1, bestJet2 = -1;
  for (int shift = 0; shift != maxshift; ++ shift) {
    bestMet = std::max(bestMet, gt->pfmet[shift]);
    bestJet1 = std::max(bestJet1, gt->jotPt[0][shift]);
    bestJet2 = std::max(bestJet2, gt->jotPt[1][shift]);
  }

  // ZnnHbb
  if (bestMet>150 && bestJet1 > 50 && bestJet2 > 25 &&
      (gt->hbbpt[0]>50 || gt->nFatjet>0)) 
  {
    return true;
  }

  // WlnHbb
  if (bestJet1>25 && bestJet2>25 &&
      ((gt->nTightElectron>0 && gt->electronPt[0]>30) ||
       (gt->nTightMuon>0 && gt->muonPt[0]>25)) &&
      (gt->hbbpt[0]>50 || gt->nFatjet>0)) 
  {
    return true;
  }

  // ZllHbb
  if (bestJet1>25 && bestJet2>25 && 
      ((gt->nLooseElectron>1 && gt->electronPt[0]>20 && gt->electronPt[1]>15) ||
       (gt->nLooseMuon>1 && gt->muonPt[0]>20 && gt->muonPt[1]>10)) &&
      (gt->hbbpt[0]>0 || gt->nFatjet>0))
  {
    return true;
  }

  return false;
}
