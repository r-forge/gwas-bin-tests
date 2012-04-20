/* This file was automatically generated using perl script ./generate_sfmt_sources.pl. Do not edit ! */

#include "config.h"

#ifdef HAVE_SSE2

#undef HAVE_SSE2
#undef HAVE_ALTIVEC

#define HAVE_SSE2
#define MEXP 19937

#include "SFMT_SSE2_19937.hpp"
namespace SFMT_SSE2_19937 {
#include "RNG.cpp.t"
}
#endif
