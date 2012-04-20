/* This file was automatically generated using perl script ./generate_sfmt_sources.pl. Do not edit ! */

#include "config.h"

#ifdef HAVE_SSE2

#undef HAVE_SSE2
#undef HAVE_ALTIVEC

#define HAVE_SSE2
#define MEXP 11213

#include "SFMT_SSE2_11213.hpp"
namespace SFMT_SSE2_11213 {
#include "RNG.cpp.t"
}
#endif
