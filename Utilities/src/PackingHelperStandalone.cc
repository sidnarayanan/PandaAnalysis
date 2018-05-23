#include "../interface/PackingHelperStandalone.h"

#include <cmath>

Double_t
pa::PackingHelperStandalone::up(UShort_t p) // unpackUnbound
{
  if(p==0) return 0;
  union {
    float flt;
    uint32_t i32;
  } conv;

  UShort_t offset;
  UInt_t mantissa,exponent;
  Int_t i;
  
  // set offset
  if(p>>10 == 0 || p>>10 == 32) offset=0;
  else                          offset=1024;
  
  // set mantissa
  i=offset+(p&0x3ff);
  unsigned m = (i<<13);
  unsigned e = 0;
  while ((m & 0x00800000) == 0) { // While not normalized
    e -= 0x00800000; // Decrement exponent (1<<23)
    m <<= 1; // Shift mantissa
  }
  m &= ~0x00800000; // Clear leading 1 bit
  e += 0x38800000; // Adjust bias ((127-14)<<23)
  mantissa = m | e; 
  
  // set exponent
  i=p>>10;
  if(i==0)
    exponent=0;
  else if(i<=30)
    exponent = (i<<23);
  else if(i==32)
    exponent = 0x47800000;
  else if(i==33)
    exponent = 0x80000000u;
  else if(i<=62)
    exponent = 0x80000000u | ((i - 32) << 23);
  else if(i==63)
    exponent = 0xC7800000;
  
  conv.i32 = mantissa+exponent;
  return conv.flt;
}

Double_t
pa::PackingHelperStandalone::up8(Char_t i, Double_t min, Double_t max, UChar_t baseminus1) // unpack8LogBound
{
  if (baseminus1 > 127)
    baseminus1 = 127;

  double l;
  if (std::abs(i) == baseminus1)
    l = max;
  else
    l = min + std::abs(i) * (max - min) / baseminus1;

  double val(std::exp(l));
  if (i < 0)
    return -val;
  else
    return val;
}
