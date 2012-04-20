/* This file was automatically generated using perl script ./generate_sfmt_sources.pl. Do not edit ! */

#include "config.h"

#ifdef HAVE_AVX

#undef HAVE_SSE2
#undef HAVE_ALTIVEC

#define HAVE_ALTIVEC
#define MEXP 44497

#include "SFMT_ALTIVEC_44497.hpp"
namespace SFMT_ALTIVEC_44497 {
#include "RNG.cpp.t"
}
#endif
