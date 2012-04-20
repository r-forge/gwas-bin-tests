/* This file was automatically generated using perl script ./generate_sfmt_sources.pl. Do not edit ! */

#include "config.h"

#ifdef ALWAYS

#undef HAVE_SSE2
#undef HAVE_ALTIVEC

#define HAVE_NORMAL
#define MEXP 19937

#include "SFMT_NORMAL_19937.hpp"
namespace SFMT_NORMAL_19937 {
#include "RNG.cpp.t"
}
#endif
