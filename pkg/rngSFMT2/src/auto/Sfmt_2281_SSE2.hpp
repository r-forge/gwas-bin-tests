/* This file was automatically generated using perl script ./generate_sfmt_classes.pl. Do not edit ! */
#ifndef SFMT_2281_SSE2_HPP
#define SFMT_2281_SSE2_HPP

#define ALWAYS
#ifdef HAVE_SSE2 // otherwise not supported

// force neutral implementation
#ifdef HAVE_SSE2
#define RESET_HAVE_SSE2
#undef HAVE_SSE2
#endif

#ifdef HAVE_ALTIVEC
#define RESET_HAVE_ALTIVEC
#undef HAVE_ALTIVEC
#endif

#define HAVE_SSE2 1

#ifndef READ_CONFIG_H
#include "../config.h"
#define READ_CONFIG_H
#endif

#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#ifdef HAVE_EMMINTRIN_H
#include <emmintrin.h>
#endif

#ifdef HAVE_ALTIVEC_H
#include <altivec.h>
#endif

#define MEXP 2281


#include <string>
#include <string.h>
#include <assert.h>
namespace SFMT_2281_SSE2 {
#include "../Sfmt_base.hpp"
}

#undef HAVE_SSE2

#ifdef RESET_HAVE_SSE2
#define HAVE_SSE2
#endif

#ifdef RESET_HAVE_ALTIVEC
#define HAVE_ALTIVEC_H
#endif

#undef SFMT_PARAMS2281_H

#undef MEXP
#undef SFMT_H
#undef SFMT_PARAMS_H
#undef SFMT_ALTI_H
#undef SFMT_SSE2_H
#undef N
#undef N32
#undef N64
#undef POS1
#undef SL1
#undef SL2
#undef SR1
#undef SR2
#undef MSK1
#undef MSK2
#undef MSK3
#undef MSK4
#undef PARITY1
#undef PARITY2
#undef PARITY3
#undef PARITY4
#undef ALTI_SL1
#undef ALTI_SR1
#undef ALTI_MSK
#undef ALTI_MSK64
#undef ALTI_SL2_PERM
#undef ALTI_SL2_PERM64
#undef ALTI_SR2_PERM
#undef ALTI_SR2_PERM64
#undef IDSTR
#endif // when
#endif // PROTECT
