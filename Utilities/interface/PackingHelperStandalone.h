#ifndef PandaTree_Objects_PackingHelperStandalone_h
#define PandaTree_Objects_PackingHelperStandalone_h

#include "Rtypes.h"

// This standalone class allows you to unpack interactively instead of relying on tabulated info.
// After compiling PandaAnalysis using scram, you can call the following in CINT:
//   #include "PandaAnalysis/Utilities/src/PackingHelperStandalone.cc"
// and then perform a TTree Draw or Scan e.g.
//   events->Scan("panda::PHS::up(genParticles.packedPt)")
// This is the analog of panda::PackingHelper::unpackUnbound.
// For the analog of panda::PackingHelper::unpack8LogBound, use panda::PHS::up8(...)


namespace pa {
  class PackingHelperStandalone {
  public:
    PackingHelperStandalone() { }

    static Double_t up(UShort_t);
    static Double_t up8(Char_t, Double_t min, Double_t max, UChar_t baseminus1);
  };
  typedef PackingHelperStandalone PHS;
}

#endif
