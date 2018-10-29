#include "../interface/CatMods.h"
#include <algorithm>

using namespace pa; 
using namespace std; 

int VBFCatMod::categorize() 
{
  if (gt.nLoosePhoton>0 || gt.nTau>0)
    return 0; 
 
  if (gt.nJot[0] < 2) 
    return 0; 

  int tight = (gt.jotEta[0]*gt.jotEta[1]<0 && 
               gt.jotPt[0][0]>80 && gt.jotPt[0][1]>40 && 
               fabs(gt.jotEta[0])<4.7 && fabs(gt.jotEta[1])<4.7) 
              ? 1 : -1; 

 float minU = min(gt.pfmet[0], min(gt.pfUWmag[0], gt.pfUZmag[0]));
  if (minU < 200) 
    return 0; 

  if (fabs(gt.calomet - minU) / minU > 0.5)
    return 0;   


  if (gt.nLooseMuon==0 && gt.nLooseElectron==0) {
    if (gt.pfmet[0] > 200)
      return 1 * tight; 
  } else if (gt.nLooseElectron==0 && gt.nTightMuon==1) {
    if (gt.nLooseMuon==1 && gt.pfUWmag[0]>200 && gt.mT[0]<160) 
      return 2 * tight;
    else if (gt.nLooseMuon==2 && gt.pfUZmag[0]>200 && 
             60 < gt.diLepMass && gt.diLepMass<120)
      return 4 * tight;
  } else if (gt.nLooseMuon==0 && gt.nTightElectron==1) {
    if (gt.nLooseElectron==1 && gt.pfUWmag[0]>200 && gt.mT[0]<160 && gt.pfmet[0]>60)
      return 3 * tight; 
    else if (gt.nLooseElectron==2 && gt.pfUZmag[0]>200 && 
             60 < gt.diLepMass && gt.diLepMass<120)
      return 5 * tight; 
  }

  return 0;
} 
