/**
 * \file EnergyCorrelations.h
 * \brief Optimized code to calculate energy correlation functions. Based on code from fj-contrib
 * \author S.Narayanan
 */
#include "fastjet/PseudoJet.hh"
#include <vector>
#include <utility>
#include <tuple>
#include <map>
#include <type_traits>
#include "boost/multiprecision/gmp.hpp"

#include "TMath.h"
#include "TString.h"
#include "PandaCore/Tools/interface/Common.h"

#ifndef PANDA_ECF_H
#define PANDA_ECF_H

// namespace mp = boost::multiprecision;
// typedef  mp::mpf_float_50 mpfloat; 
typedef double mpfloat;
template <typename D>
inline mpfloat d2m(D x) { return static_cast<mpfloat>(x); }


namespace pandaecf {
  /**
   * \brief delta-r-squared metric between two pseudojets
   * @param  j1 first jet
   * @param  j2 second jet
   * @return    \f$dR^2\f$
   */
   mpfloat DeltaR2(const fastjet::PseudoJet& j1, const fastjet::PseudoJet& j2);

  class Calculator {
  public:
    enum param {
      oP=0, nP, bP, ecfP
    };
    typedef std::tuple<int, int, int, mpfloat> data_type;
    typedef std::tuple<int, int, int> pos_type;

    Calculator(int maxN = 4,
               std::vector<float> bs = {0.5, 1, 2, 4});
    ~Calculator() { }

    data_type access(int pos) const { return access(_oneToThree(pos)); }
    data_type access(pos_type pos) const;
    void calculate(const std::vector<fastjet::PseudoJet>&);

    // just a forward iterator
    class iterator {
    private:
      const Calculator *_c;
      int _pos;
      data_type _data;

      void _access() { _data = _c->access(_pos); }
    public:
      iterator(const Calculator *c, int pos = 0): _c(c), _pos(pos) { _access(); }
      iterator(const iterator& rhs): _c(rhs._c), _pos(rhs._pos), _data(rhs._data) { }
      ~iterator() { }

      iterator& operator++() { _pos++; _access(); return *this; }
      iterator operator++(int) { auto old(*this); ++(*this); return old; }
      iterator operator+(int n) const { return iterator(_c, _pos+n); }
      iterator& operator+=(int n) { _pos += n; _access(); return *this; }
      int operator-(const iterator& rhs) const { return this->_pos - rhs._pos; }
      bool operator==(const iterator& rhs) const { return this->_pos == rhs._pos; }
      bool operator!=(const iterator& rhs) const { return !( (*this) == rhs ); }
      const data_type& operator->() const { return _data; }
      template <int I>
        const typename std::tuple_element<I, data_type>::type& get() const { return std::get<I>(_data); }
        // leaving this here just because I'm kind of proud of this C++11 mess
        // auto get() -> typename std::add_const<typename std::add_lvalue_reference<decltype(std::get<I> (_data))>::type>::type const { return std::get<I>(_data); }

    };

    iterator begin() const { return iterator(this, 0); }
    iterator end() const { return iterator(this, _bN * _nN * _oN); } 

  private:
    void _set(pos_type pos, mpfloat x) { _ecfs[_threeToOne(pos)] = x; }
    pos_type _oneToThree(int pos) const;
    int _threeToOne(pos_type pos) const { return std::get<oP>(pos) 
                                                 + _oN * std::get<nP>(pos) 
                                                 + _oN + _nN * std::get<bP>(pos); }

    std::vector<float> _bs;
    std::vector<int> _ns, _os;
    const int _bN, _nN, _oN;
    std::vector<mpfloat> _ecfs;
    std::vector<mpfloat> pT; // these are member variables just to avoid
    std::vector<std::vector<mpfloat>> dR, dRBeta; // re-allocating memory 
  };
}

#endif
