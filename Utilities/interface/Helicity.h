#include <TLorentzVector.h>
#include <TVector3.h>
// Transform to Collins-Soper frame then compute the
// angle between the lepton momentum and an axis
// that bisects the angle between the quark and opposite
// to the anti quark direction
double CosThetaCollinsSoper(TLorentzVector LVlep1, TLorentzVector LVlep2) {
  TLorentzVector LVZ = LVlep1+LVlep2;
  //do transformation to Collins-Soper frame (1 rotation, 2 boosts)
  // 1st transormation - rotate to ptZ direction
  double zrot = -LVZ.Phi();
  LVlep1.RotateZ(zrot);
  LVlep2.RotateZ(zrot);
  
  // 2nd transformation - boost in z direction
  double beta_boostz = -LVZ.Pz()/LVZ.E();
  LVlep1.Boost(0.,0.,beta_boostz);
  LVlep2.Boost(0.,0.,beta_boostz);
  
  // 3rd transformation: boost in transverse direction (x-prime)
  double beta_boostx = -(LVlep1.Px()+LVlep2.Px())/(LVlep1.E()+LVlep2.E());
  LVlep1.Boost(beta_boostx,0.,0.);
  LVlep2.Boost(beta_boostx,0.,0.);
  
  // compute cos(theta*) in Colin-Soper frame
  double cos_theta_CS = LVlep1.CosTheta();
  if (LVZ.Pz() < 0) {
         cos_theta_CS *= -1.;
  }
  
  if     (cos_theta_CS >= +1) cos_theta_CS = +0.999999999;
  else if(cos_theta_CS <= -1) cos_theta_CS = -0.999999999;
  
  return cos_theta_CS;
}

// Helicity angles
double CosThetaStar(TLorentzVector lep1, TLorentzVector lep2, TLorentzVector grandma) {
  // Leptons in laboratory frame
  TLorentzVector p4_l1_mother = lep1;
  TLorentzVector dilep=lep1+lep2;
  // Boost back to dilepton (mother) rest frame
  p4_l1_mother.Boost(-dilep.X()/dilep.T(), -dilep.Y()/dilep.T(), -dilep.Z()/dilep.T());
  // Z in laboratory frame
  TLorentzVector p4_mother_grandma = dilep;
  // Boost back to Z+MET (grandmother) rest frame
  p4_mother_grandma.Boost( -grandma.X()/grandma.T(), -grandma.Y()/grandma.T(), -grandma.Z()/grandma.T());
  // Get the 3 vectors of the lepton in the mother's rest frame, and the mother in the grandmother's rest frame
  TVector3 p3_mother_grandma = p4_mother_grandma.Vect(), p3_l1_mother = p4_l1_mother.Vect();

  double cos_theta_star_l1 = p3_mother_grandma.Dot(p3_l1_mother) / ( p3_mother_grandma.Mag() * p3_l1_mother.Mag() );
  return cos_theta_star_l1;
}
