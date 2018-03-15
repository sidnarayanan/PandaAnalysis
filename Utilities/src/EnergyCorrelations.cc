/**
 * \file EnergyCorrelations.cc
 * \brief Optimized code to calculate energy correlation functions. Based on code from fj-contrib
 * \author S.Narayanan
 */
#include "../interface/EnergyCorrelations.h"
#define PI 3.141592654

using namespace pandaecf;
using namespace std;
typedef Calculator C;

mpfloat pandaecf::DeltaR2(const fastjet::PseudoJet& j1, const fastjet::PseudoJet& j2) 
{
    mpfloat dEta{j1.eta()-j2.eta()}; 
    mpfloat dPhi{j1.phi()-j2.phi()};

    if (dPhi<-PI)
        dPhi = 2*PI+dPhi;
    else if (dPhi>PI)
        dPhi = -2*PI+dPhi;

    return dEta*dEta + dPhi*dPhi;
}

C::Calculator(int maxN, vector<float> bs):
  _bs(bs),
  _os({1,2,3}),
  _bN(_bs.size()),
  _nN(maxN),
  _oN(_os.size())
{
  _ns = vector<int>(maxN);
  for (int i = 0; i != maxN; ++i)
    _ns[i] = i + 1;
  _ecfs.resize(_nN * _oN * _bN);
}

C::data_type C::access(C::pos_type pos) const
{
  if (_threeToOne(pos) == _oN * _nN * _bN) {
    // don't throw error - fail silently, to allow end() to be defined
    return make_tuple(-1,-1,-1,d2m(-1));
  }
  return make_tuple(get<oP>(pos),
                    get<nP>(pos),
                    get<bP>(pos),
                    _ecfs.at(_threeToOne(pos)));
}

C::pos_type C::_oneToThree(int pos) const
{
  auto oI = pos % _oN;
  pos -= oI; pos /= _oN;
  auto nI = pos % _nN;
  pos -= nI; pos /= _nN;
  auto bI = pos;
  return make_tuple(oI, nI, bI);
}

void C::calculate(const vector<fastjet::PseudoJet> &particles)
{
  if (_nN == 0 || _bN == 0 || _oN == 0)
    return;

  int nParticles = particles.size();

  // cache kinematics
  pT.resize(nParticles, d2m(0)); dR.resize(nParticles); dRBeta.resize(nParticles);
  for (int iP=0; iP!=nParticles; ++iP) {
    dR[iP].resize(nParticles, d2m(0));
    dRBeta[iP].resize(nParticles, d2m(0));
  }

  for (int iP=0; iP!=nParticles; ++iP) {
    const fastjet::PseudoJet& pi = particles[iP];
    pT[iP] = pi.perp();
    for (int jP=0; jP!=iP; ++jP) {
      if (iP == jP)
        continue;
      const fastjet::PseudoJet& pj = particles[jP];
      dR[iP][jP] = DeltaR2(pi,pj);
    }
  }

  for (int bI = 0; bI != _bN; ++bI) {
    // get the normalization factor
    mpfloat baseNorm{0};
    for (auto& pt : pT)
      baseNorm += pt;
    
    // reweight angles
    mpfloat halfBeta{_bs.at(bI) / 2.};
    for (int iP=0; iP!=nParticles; ++iP) {
      for (int jP=0; jP!=iP; ++jP) {
        if (iP == jP)
          continue;
        dRBeta[iP][jP] = pow(dR[iP][jP], halfBeta);
      }
    }

    // now we compute the ECFNs
    // trivial case, n = 1
    for (int oI = 0; oI != _oN; ++oI) {
      _set(make_tuple(oI, 0, bI), d2m(1));
    }

    if (_nN < 2)
      return;

    mpfloat norm2{pow(baseNorm, 2)};
    mpfloat norm3{pow(baseNorm, 3)};
    mpfloat norm4{pow(baseNorm, 4)};
    vector<vector<mpfloat>> vals(_nN);
    for (int nI = 0; nI != _nN; ++nI) {
      vals[nI].resize(_oN, d2m(0));
    }

    // now we loop
    vector<mpfloat> angles3(3), angles4(6); // N(N-1)/2
    for (int iP = 0; iP != nParticles; ++iP) {
      for (int jP = 0; jP != nParticles; ++jP) {
        const mpfloat pt_ij{pT[iP] * pT[jP]};
        const mpfloat angle_ij{dRBeta[iP][jP]};

        angles3[0] = angle_ij; angles4[0] = angle_ij;
        vals[1][0] += pt_ij * angle_ij;

        if (_nN > 2) { 
          for (int kP = 0; kP != nParticles; ++kP) {
            const mpfloat angle_ik{dRBeta[iP][kP]}, angle_jk{dRBeta[jP][kP]};
            const mpfloat pt_ijk{pt_ij * pT[kP]};

            angles3[1] = angle_ik; angles3[2] = angle_jk;
            angles4[1] = angle_ik; angles4[2] = angle_jk;

            insertion_sort(angles3);
            
            mpfloat inc{pt_ijk};
            for (int oI = 0; oI != _oN; ++oI) {
              inc *= angles3[oI];
              vals[2][oI] += inc;
            }

            if (_nN > 3) {
              for (int lP = 0; lP != nParticles; ++lP) {
                const mpfloat pt_ijkl{pt_ijk * pT[lP]};
                angles4[3] = dRBeta[iP][lP]; 
                angles4[4] = dRBeta[jP][lP];
                angles4[5] = dRBeta[kP][lP];

                partial_sort(angles4.begin(), angles4.begin() + 3, angles4.end());

                mpfloat inc{pt_ijkl};
                for (int oI = 0; oI != _oN; ++oI) {
                  inc *= angles4[oI];
                  vals[3][oI] += inc;
                }
                
              } // l
            } // N > 3

          } // k
        } // N > 2
      } // j
    } // i

  
    // set the values
    mpfloat val2 = vals[1][0];
    val2 /= norm2;
    for (int oI = 0; oI != _oN; ++oI) {
      _set(make_tuple(oI, 1, bI), val2);
      _set(make_tuple(oI, 2, bI), vals[2][oI] / norm3);
      _set(make_tuple(oI, 3, bI), vals[3][oI] / norm4);
    }
  } // beta loop
}

