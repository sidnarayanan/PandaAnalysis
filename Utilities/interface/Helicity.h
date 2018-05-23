#include <TLorentzVector.h>
#include <TVector3.h>
namespace pa {
    // Transform to Collins-Soper frame then compute the
    // angle between the lepton momentum and an axis
    // that bisects the angle between the quark and opposite
    // to the anti quark direction
    double CosThetaCollinsSoper(TLorentzVector LVlep1, TLorentzVector LVlep2); 

    // Helicity angles
    double CosThetaStar(TLorentzVector lep1, TLorentzVector lep2, TLorentzVector grandma); 
}
