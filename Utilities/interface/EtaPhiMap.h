#ifndef CROMBIE_ETAPHIMAP_H
#define CROMBIE_ETAPHIMAP_H 1

#include <vector>
#include <cmath>
#include <functional>
#include "TMath.h"
#include "TVector2.h"
#include "PandaCore/Tools/interface/Common.h"
namespace pa {
  template<typename T>
  class EtaPhiMap {
   public:
    EtaPhiMap (double spacing, double maxeta = 5.0,
               std::function<double(const T*)> eta = [](const T* p){ return p->eta(); },
               std::function<double(const T*)> phi = [](const T* p){ return p->phi(); })
      : _spacing{spacing}, _maxeta{maxeta},
        n_etabins{2 * (static_cast<unsigned>(std::ceil(_maxeta/_spacing)) + 1)}, // Add one for over/underflow and times two for positive and negative
        n_phibins{static_cast<unsigned>(std::ceil(TMath::TwoPi()/_spacing))},    // Always transform phi into [0, 2pi) in bin
        geteta{eta}, getphi{phi} {
  
          particles.resize(n_etabins * n_phibins);
  
        }
  
    /// Add a collection of particles to the grid. Call once per event, because this gets cleared
    template<typename C> void AddParticles (C& collection);
  
    /// Get particles within dr of a given eta, phi
    std::vector<const T*> GetParticles(double eta, double phi, double dr);
  
   private:
  
    /// Stores all of the particles
    std::vector<std::vector<const T*>> particles;
    /// Spacing of each grid
    double _spacing;
    /// Max eta of grid (has overflow and underflow bins)
    double _maxeta;
  
    /// Number of eta bins
    unsigned n_etabins;
    /// Number of phi bins
    unsigned n_phibins;
  
    std::function<double(const T*)> geteta;  ///< Function for getting eta from pointers
    std::function<double(const T*)> getphi;  ///< Function for getting phi from pointers
  
    /// Reset the particles in each grid point
    void clear();
    /// Get the bin number
    unsigned bin(double eta, double phi);
    unsigned bin(unsigned eta_bin, unsigned phi_bin);
    /// Get eta bin
    unsigned etabin(double eta);
    /// Get phi bin
    unsigned phibin(double phi);
  
  };
  
  
  template<typename T>
  template<typename C>
  void EtaPhiMap<T>::AddParticles (C& collection) {
  
    // Only call this function when adding a full collection
    clear();
  
    for (auto& cand : collection) {
      auto* ptr = &cand;
      particles[bin(geteta(ptr), getphi(ptr))].push_back(ptr);
    }
  
  }
  
  
  template<typename T>
  std::vector<const T*> EtaPhiMap<T>::GetParticles(double eta, double phi, double dr) {
    double dr2 = std::pow(dr, 2);
  
    std::vector<const T*> output;
  
    auto min_eta = etabin(eta - dr);
    auto max_eta = etabin(eta + dr);
    auto min_phi = phibin(phi - dr);
    auto max_phi = phibin(phi + dr);
    auto running_phi = min_phi;
    while (true) {
      for (auto running_eta = min_eta; running_eta <= max_eta; ++running_eta) {
        for (auto* particle : particles[bin(running_eta, running_phi)]) {
          if (DeltaR2(eta, phi, geteta(particle), getphi(particle)) < dr2)
            output.push_back(particle);
        }
      }
      // If at the end of phi, break after processing
      if (running_phi == max_phi)
        break;
      // Increment and go to zero if at the end of phi
      if (++running_phi == n_phibins)
        running_phi = 0;
    }
  
    return output;
  
  }
  

  template<typename T>
  unsigned EtaPhiMap<T>::bin(double eta, double phi) {
    return bin(etabin(eta), phibin(phi));
  }
  
  template<typename T>
  unsigned EtaPhiMap<T>::bin(unsigned eta_bin, unsigned phi_bin) {
    return n_etabins * phi_bin + eta_bin;
  }
  
  template<typename T>
  unsigned EtaPhiMap<T>::etabin(double eta) {
    return std::max(0u, std::min(n_etabins - 1, n_etabins/2 + static_cast<unsigned>(eta/_spacing)));
  }
  
  template<typename T>
  unsigned EtaPhiMap<T>::phibin(double phi) {
    return TVector2::Phi_0_2pi(phi)/_spacing;
  }
  
  
  template<typename T>
  void EtaPhiMap<T>::clear() {
    for (auto& grid : particles)
      grid.clear();
  }
  
}

#endif

