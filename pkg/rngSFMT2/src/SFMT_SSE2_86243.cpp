/* This file was automatically generated using perl script ./generate_sfmt_sources.pl. Do not edit ! */

#include "config.h"

#ifdef HAVE_SSE2

#undef HAVE_SSE2
#undef HAVE_ALTIVEC

#define HAVE_SSE2
#define MEXP 86243

#include "SFMT_SSE2_86243.hpp"
namespace SFMT_SSE2_86243 {
#include "RNG.cpp.t"
}
#endif
