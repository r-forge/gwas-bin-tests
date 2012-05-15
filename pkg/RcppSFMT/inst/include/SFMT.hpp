#ifndef SFMT_HPP
#define SFMT_HPP

// CRITICAL: set MEXP and the SIMD implementation type (HAVE_SSE2, HAVE_ALTIVEC), etc...
#include "config.h"

// The configure option --disable-sse2 set HAVE_SSE2 to 0 but does not undef it
// because I was not able to dot it. So here we check this value
#ifdef HAVE_SSE2
#if HAVE_SSE2 == 0
#undef HAVE_SSE2
#warning undefining HAVE_SSE2
#endif
#endif

// same for --disable-altivec and HAVE_ALTIVEC
#ifdef HAVE_ALTIVEC
#if HAVE_ALTIVEC == 0
#undef HAVE_ALTIVEC
#warning undefining HAVE_ALTIVEC
#endif
#endif

// include the defined constants of the original C implementations as well as some functions
extern "C" {
#define static // trick to remove the C static stuff
#include "SFMT.c"
#undef static
}

#include <string>
#include <vector>
namespace RNG {
using namespace std;

class SFMT {
private: // ========== DECLARATION OF INSTANCE STATE VARIABLES =========

		// ===== instance variable ====
		// WRAP INSTRUCTIONS 1:
		// These state variables come from the "FILE GLOBAL VARIABLES" section
		// of SFMT.c: take all variables without the "static" keyword
		// and put the definitions in the init() method
		// because it is not allowed  in C++ to define vars inside a class


		// ========== BEGIN ===========

		/** the 128-bit internal state array */
		w128_t sfmt[N];
		/** the 32bit integer pointer to the 128-bit internal state array */
		uint32_t *psfmt32; // = &sfmt[0].u[0];
		#if !defined(BIG_ENDIAN64) || defined(ONLY64)
		/** the 64bit integer pointer to the 128-bit internal state array */
		uint64_t *psfmt64; // = (uint64_t *)&sfmt[0].u[0];
		#endif
		/** index counter to the 32-bit internal state array */
		int idx;
		/** a flag: it is 0 if and only if the internal state is not yet
		 * initialized. */
		int initialized; // = 0;
		/** a parity check vector which certificate the period of 2^{MEXP} */
		uint32_t parity[4]; // = {PARITY1, PARITY2, PARITY3, PARITY4};

		// ========== END ===========

	public: // ========== LIFECYCLE ======================
		/**
		 * init the rng using SFMT C function gen_rand_all() :
		 * This function fills the internal state array with pseudorandom
		 * integers.
		 */
		SFMT() { init(); gen_rand_all(); }

		/**
		 * init the rng with a seed
		 * @see set_seed
		 */
		SFMT(uint32_t seed) { init(); set_seed(seed); }

		/**
		 * init the rng with a two seeds
		 * @see set_seed
		 */
		SFMT(uint32_t seed1, uint32_t seed2) {
			init();
			set_seeds(seed1, seed2);
		}

		/**
		 * init the rng with an array of seeds
		 * @see set_seeds
		 */
		SFMT(uint32_t *init_key, int key_length) {
			init();
			set_seeds(init_key, key_length);
		}



	public: // ===== PUBLIC INTERFACE =====
		/**
		 * This function generates and returns 32-bit pseudorandom number.
		 *
		 */
		uint32_t rand32() { return gen_rand32(); }

		/**
		 *
		 */
		void rand32(vector<uint32_t>& v) {
			int nb = v.size();
			if ( nb == 0 )
				return;
			fill_array32(&v[0], nb);
		}


		/**
		 * generates a random number on (0,1)-real-interval
		 *
		 * N.B: call genrand_real3()
		 */
		double random() { return genrand_real3(); }

		/**
		 * set the seed. Can be used at any time to reset the seed
		 *
		 * Call the init_gen_rand() SFMT function:
		 *
		 * This function initializes the internal state array with a 32-bit
		 * integer seed.
		 *
		 * @param seed a 32-bit integer used as the seed.
		 */
		void set_seed(uint32_t seed) { init_gen_rand(seed); }

		/**
		 * set the seeds. Can be used at any time to reset the seed
		 *
		 * Just a wrapper around the init_by_array() SFMT funbction:
		 *
		 * This function initializes the internal state array,
		 * with two seeds
		 *
		 * @param seed1 a 32-bit integer used as the first seed.
		 * @param seed2 a 32-bit integer used as the ssecond eed.
		 */
		 void set_seeds(uint32_t seed1, uint32_t seed2) {
				uint32_t seeds[] = { seed1, seed2 };
				set_seeds(seeds, 2);
		 }

		/**
		 * set the seeds. Can be used at any time to reset the seed
		 *
		 * Just a wrapper around the init_by_array() SFMT funbction:
		 *
		 * This function initializes the internal state array,
		 * with an array of 32-bit integers used as the seeds
		 * @param init_key the array of 32-bit integers, used as a seed.
		 * @param key_length the length of init_key.
		 */
		 void set_seeds(uint32_t *init_key, int key_length) { init_by_array(init_key, key_length);  }

		 /**
		  * get the Mersenne Exponent used by this implementation
		  *
		  * @return the exponent, as an integer
		  */
		 static int GET_MERSENNE_EXPONENT() {
			 return MEXP;
		 }


		static std::string ID()  {
			return get_idstring();
		}

	protected:
		/**
		 * Init state variables: initialization as copy/paste
		 * of section "FILE GLOBAL VARIABLES" in SFMT.c
		 *
		 */
		 void init() {

			// WRAP INSTRUCTIONS 2:
			// This corresponds to the initialization of state variables that comes
			 // from the "FILE GLOBAL VARIABLES" section  of SFMT.c

			// ========== BEGIN ===========
			/** the 32bit integer pointer to the 128-bit internal state array */
			psfmt32 = &sfmt[0].u[0];
			#if !defined(BIG_ENDIAN64) || defined(ONLY64)
			/** the 64bit integer pointer to the 128-bit internal state array */
			psfmt64 = (uint64_t *)&sfmt[0].u[0];
			#endif
			/** a flag: it is 0 if and only if the internal state is not yet
			 * initialized. */
			initialized = 0;
			/** a parity check vector which certificate the period of 2^{MEXP} */
		    parity[0] = PARITY1;
		    parity[1] = PARITY2;
		    parity[2] = PARITY3;
		    parity[3] = PARITY4;

		    // ========== END ===========

		}

	private: // ==== original inline static C SFMT functions =====

// little trick: remove the static statements
#define static

#ifdef HAVE_SSE2
// WRAP INSTRUCTIONS 3:
// copy/paste here the definitions of functions from SFMT-sse2.h: mm_recursion, 	 gen_rand_all, gen_rand_array


		 /**
		  * This function represents the recursion formula.
		  * @param a a 128-bit part of the interal state array
		  * @param b a 128-bit part of the interal state array
		  * @param c a 128-bit part of the interal state array
		  * @param d a 128-bit part of the interal state array
		  * @param mask 128-bit mask
		  * @return output
		  */
		 PRE_ALWAYS static __m128i mm_recursion(__m128i *a, __m128i *b,
		 				   __m128i c, __m128i d, __m128i mask) {
		     __m128i v, x, y, z;

		     x = _mm_load_si128(a);
		     y = _mm_srli_epi32(*b, SR1);
		     z = _mm_srli_si128(c, SR2);
		     v = _mm_slli_epi32(d, SL1);
		     z = _mm_xor_si128(z, x);
		     z = _mm_xor_si128(z, v);
		     x = _mm_slli_si128(x, SL2);
		     y = _mm_and_si128(y, mask);
		     z = _mm_xor_si128(z, x);
		     z = _mm_xor_si128(z, y);
		     return z;
		 }

		 /**
		  * This function fills the internal state array with pseudorandom
		  * integers.
		  */
		 inline static void gen_rand_all(void) {
		     int i;
		     __m128i r, r1, r2, mask;
		     mask = _mm_set_epi32(MSK4, MSK3, MSK2, MSK1);

		     r1 = _mm_load_si128(&sfmt[N - 2].si);
		     r2 = _mm_load_si128(&sfmt[N - 1].si);
		     for (i = 0; i < N - POS1; i++) {
		 	r = mm_recursion(&sfmt[i].si, &sfmt[i + POS1].si, r1, r2, mask);
		 	_mm_store_si128(&sfmt[i].si, r);
		 	r1 = r2;
		 	r2 = r;
		     }
		     for (; i < N; i++) {
		 	r = mm_recursion(&sfmt[i].si, &sfmt[i + POS1 - N].si, r1, r2, mask);
		 	_mm_store_si128(&sfmt[i].si, r);
		 	r1 = r2;
		 	r2 = r;
		     }
		 }

		 /**
		  * This function fills the user-specified array with pseudorandom
		  * integers.
		  *
		  * @param array an 128-bit array to be filled by pseudorandom numbers.
		  * @param size number of 128-bit pesudorandom numbers to be generated.
		  */
		 inline static void gen_rand_array(w128_t *array, int size) {
		     int i, j;
		     __m128i r, r1, r2, mask;
		     mask = _mm_set_epi32(MSK4, MSK3, MSK2, MSK1);

		     r1 = _mm_load_si128(&sfmt[N - 2].si);
		     r2 = _mm_load_si128(&sfmt[N - 1].si);
		     for (i = 0; i < N - POS1; i++) {
		 	r = mm_recursion(&sfmt[i].si, &sfmt[i + POS1].si, r1, r2, mask);
		 	_mm_store_si128(&array[i].si, r);
		 	r1 = r2;
		 	r2 = r;
		     }
		     for (; i < N; i++) {
		 	r = mm_recursion(&sfmt[i].si, &array[i + POS1 - N].si, r1, r2, mask);
		 	_mm_store_si128(&array[i].si, r);
		 	r1 = r2;
		 	r2 = r;
		     }
		     /* main loop */
		     for (; i < size - N; i++) {
		 	r = mm_recursion(&array[i - N].si, &array[i + POS1 - N].si, r1, r2,
		 			 mask);
		 	_mm_store_si128(&array[i].si, r);
		 	r1 = r2;
		 	r2 = r;
		     }
		     for (j = 0; j < 2 * N - size; j++) {
		 	r = _mm_load_si128(&array[j + size - N].si);
		 	_mm_store_si128(&sfmt[j].si, r);
		     }
		     for (; i < size; i++) {
		 	r = mm_recursion(&array[i - N].si, &array[i + POS1 - N].si, r1, r2,
		 			 mask);
		 	_mm_store_si128(&array[i].si, r);
		 	_mm_store_si128(&sfmt[j++].si, r);
		 	r1 = r2;
		 	r2 = r;
		     }
		 }


#endif


#ifdef HAVE_ALTIVEC
// WRAP INSTRUCTIONS 4:
// copy/paste here the definitions of functions from SFMT-alti.h:  gen_rand_all, gen_rand_array, swap

/**
* This function represents the recursion formula in AltiVec and BIG ENDIAN.
* @param a a 128-bit part of the interal state array
* @param b a 128-bit part of the interal state array
* @param c a 128-bit part of the interal state array
* @param d a 128-bit part of the interal state array
* @return output
*/
inline static vector unsigned int vec_recursion(vector unsigned int a,
					vector unsigned int b,
					vector unsigned int c,
					vector unsigned int d) {

 const vector unsigned int sl1 = ALTI_SL1;
 const vector unsigned int sr1 = ALTI_SR1;
#ifdef ONLY64
 const vector unsigned int mask = ALTI_MSK64;
 const vector unsigned char perm_sl = ALTI_SL2_PERM64;
 const vector unsigned char perm_sr = ALTI_SR2_PERM64;
#else
 const vector unsigned int mask = ALTI_MSK;
 const vector unsigned char perm_sl = ALTI_SL2_PERM;
 const vector unsigned char perm_sr = ALTI_SR2_PERM;
#endif
 vector unsigned int v, w, x, y, z;
 x = vec_perm(a, (vector unsigned int)perm_sl, perm_sl);
 v = a;
 y = vec_sr(b, sr1);
 z = vec_perm(c, (vector unsigned int)perm_sr, perm_sr);
 w = vec_sl(d, sl1);
 z = vec_xor(z, w);
 y = vec_and(y, mask);
 v = vec_xor(v, x);
 z = vec_xor(z, y);
 z = vec_xor(z, v);
 return z;
}

/**
* This function fills the internal state array with pseudorandom
* integers.
*/
inline static void gen_rand_all(void) {
 int i;
 vector unsigned int r, r1, r2;

 r1 = sfmt[N - 2].s;
 r2 = sfmt[N - 1].s;
 for (i = 0; i < N - POS1; i++) {
r = vec_recursion(sfmt[i].s, sfmt[i + POS1].s, r1, r2);
sfmt[i].s = r;
r1 = r2;
r2 = r;
 }
 for (; i < N; i++) {
r = vec_recursion(sfmt[i].s, sfmt[i + POS1 - N].s, r1, r2);
sfmt[i].s = r;
r1 = r2;
r2 = r;
 }
}

/**
* This function fills the user-specified array with pseudorandom
* integers.
*
* @param array an 128-bit array to be filled by pseudorandom numbers.
* @param size number of 128-bit pesudorandom numbers to be generated.
*/
inline static void gen_rand_array(w128_t *array, int size) {
 int i, j;
 vector unsigned int r, r1, r2;

 r1 = sfmt[N - 2].s;
 r2 = sfmt[N - 1].s;
 for (i = 0; i < N - POS1; i++) {
r = vec_recursion(sfmt[i].s, sfmt[i + POS1].s, r1, r2);
array[i].s = r;
r1 = r2;
r2 = r;
 }
 for (; i < N; i++) {
r = vec_recursion(sfmt[i].s, array[i + POS1 - N].s, r1, r2);
array[i].s = r;
r1 = r2;
r2 = r;
 }
 /* main loop */
 for (; i < size - N; i++) {
r = vec_recursion(array[i - N].s, array[i + POS1 - N].s, r1, r2);
array[i].s = r;
r1 = r2;
r2 = r;
 }
 for (j = 0; j < 2 * N - size; j++) {
sfmt[j].s = array[j + size - N].s;
 }
 for (; i < size; i++) {
r = vec_recursion(array[i - N].s, array[i + POS1 - N].s, r1, r2);
array[i].s = r;
sfmt[j++].s = r;
r1 = r2;
r2 = r;
 }
}

 #ifndef ONLY64
 #if defined(__APPLE__)
 #define ALTI_SWAP (vector unsigned char) \
	(4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11)
 #else
 #define ALTI_SWAP {4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11}
 #endif
 /**
  * This function swaps high and low 32-bit of 64-bit integers in user
  * specified array.
  *
  * @param array an 128-bit array to be swaped.
  * @param size size of 128-bit array.
  */
 inline static void swap(w128_t *array, int size) {
	 int i;
	 const vector unsigned char perm = ALTI_SWAP;

	 for (i = 0; i < size; i++) {
	array[i].s = vec_perm(array[i].s, (vector unsigned int)perm, perm);
	 }
 }
 #endif

#endif


 // WRAP INSTRUCTIONS 5:
 // copy/paste some required functions from SFMT.h: genrand_real3 and to_real3
 /** generates a random number on (0,1)-real-interval */
 inline static double to_real3(uint32_t v)
 {
     return (((double)v) + 0.5)*(1.0/4294967296.0);
     /* divided by 2^32 */
 }

 /** generates a random number on (0,1)-real-interval */
 inline static double genrand_real3(void)
 {
     return to_real3(gen_rand32());
 }

// WRAP INSTRUCTIONS 6:
// these are the function that use state global variables
// because we have instance variables with the same name
// we can use the exact same code
// so there functions are just copy/pasted here
//
// Functions: gen_rand32, gen_rand64, fill_array32, fill_array64
// init_gen_rand, init_by_array, gen_rand_all, gen_rand_array
// So just copy/paste the function defintions from SFMT.c

#if (!defined(HAVE_ALTIVEC)) && (!defined(HAVE_SSE2))
/**
 * This function fills the internal state array with pseudorandom
 * integers.
 */
inline static void gen_rand_all(void) {
    int i;
    w128_t *r1, *r2;

    r1 = &sfmt[N - 2];
    r2 = &sfmt[N - 1];
    for (i = 0; i < N - POS1; i++) {
	do_recursion(&sfmt[i], &sfmt[i], &sfmt[i + POS1], r1, r2);
	r1 = r2;
	r2 = &sfmt[i];
    }
    for (; i < N; i++) {
	do_recursion(&sfmt[i], &sfmt[i], &sfmt[i + POS1 - N], r1, r2);
	r1 = r2;
	r2 = &sfmt[i];
    }
}

/**
 * This function fills the user-specified array with pseudorandom
 * integers.
 *
 * @param array an 128-bit array to be filled by pseudorandom numbers.
 * @param size number of 128-bit pseudorandom numbers to be generated.
 */
inline static void gen_rand_array(w128_t *array, int size) {
    int i, j;
    w128_t *r1, *r2;

    r1 = &sfmt[N - 2];
    r2 = &sfmt[N - 1];
    for (i = 0; i < N - POS1; i++) {
	do_recursion(&array[i], &sfmt[i], &sfmt[i + POS1], r1, r2);
	r1 = r2;
	r2 = &array[i];
    }
    for (; i < N; i++) {
	do_recursion(&array[i], &sfmt[i], &array[i + POS1 - N], r1, r2);
	r1 = r2;
	r2 = &array[i];
    }
    for (; i < size - N; i++) {
	do_recursion(&array[i], &array[i - N], &array[i + POS1 - N], r1, r2);
	r1 = r2;
	r2 = &array[i];
    }
    for (j = 0; j < 2 * N - size; j++) {
	sfmt[j] = array[j + size - N];
    }
    for (; i < size; i++, j++) {
	do_recursion(&array[i], &array[i - N], &array[i + POS1 - N], r1, r2);
	r1 = r2;
	r2 = &array[i];
	sfmt[j] = array[i];
    }
}
#endif

#ifndef ONLY64
		/**
		 * This function generates and returns 32-bit pseudorandom number.
		 * init_gen_rand or init_by_array must be called before this function.
		 * @return 32-bit pseudorandom number
		 */
uint32_t gen_rand32(void) {
	uint32_t r;

	assert(initialized);
	if (idx >= N32) {
	gen_rand_all();
	idx = 0;
	}
	r = psfmt32[idx++];
	return r;
}
#endif
/**
 * This function generates and returns 64-bit pseudorandom number.
 * init_gen_rand or init_by_array must be called before this function.
 * The function gen_rand64 should not be called after gen_rand32,
 * unless an initialization is again executed.
 * @return 64-bit pseudorandom number
 */
uint64_t gen_rand64(void) {
#if defined(BIG_ENDIAN64) && !defined(ONLY64)
	uint32_t r1, r2;
#else
	uint64_t r;
#endif

	assert(initialized);
	assert(idx % 2 == 0);

	if (idx >= N32) {
	gen_rand_all();
	idx = 0;
	}
#if defined(BIG_ENDIAN64) && !defined(ONLY64)
	r1 = psfmt32[idx];
	r2 = psfmt32[idx + 1];
	idx += 2;
	return ((uint64_t)r2 << 32) | r1;
#else
	r = psfmt64[idx / 2];
	idx += 2;
	return r;
#endif
}

#ifndef ONLY64
/**
 * This function generates pseudorandom 32-bit integers in the
 * specified array[] by one call. The number of pseudorandom integers
 * is specified by the argument size, which must be at least 624 and a
 * multiple of four.  The generation by this function is much faster
 * than the following gen_rand function.
 *
 * For initialization, init_gen_rand or init_by_array must be called
 * before the first call of this function. This function can not be
 * used after calling gen_rand function, without initialization.
 *
 * @param array an array where pseudorandom 32-bit integers are filled
 * by this function.  The pointer to the array must be \b "aligned"
 * (namely, must be a multiple of 16) in the SIMD version, since it
 * refers to the address of a 128-bit integer.  In the standard C
 * version, the pointer is arbitrary.
 *
 * @param size the number of 32-bit pseudorandom integers to be
 * generated.  size must be a multiple of 4, and greater than or equal
 * to (MEXP / 128 + 1) * 4.
 *
 * @note \b memalign or \b posix_memalign is available to get aligned
 * memory. Mac OSX doesn't have these functions, but \b malloc of OSX
 * returns the pointer to the aligned memory block.
 */
void fill_array32(uint32_t *array, int size) {
    assert(initialized);
    assert(idx == N32);
    assert(size % 4 == 0);
    assert(size >= N32);

    gen_rand_array((w128_t *)array, size / 4);
    idx = N32;
}
#endif

/**
 * This function generates pseudorandom 64-bit integers in the
 * specified array[] by one call. The number of pseudorandom integers
 * is specified by the argument size, which must be at least 312 and a
 * multiple of two.  The generation by this function is much faster
 * than the following gen_rand function.
 *
 * For initialization, init_gen_rand or init_by_array must be called
 * before the first call of this function. This function can not be
 * used after calling gen_rand function, without initialization.
 *
 * @param array an array where pseudorandom 64-bit integers are filled
 * by this function.  The pointer to the array must be "aligned"
 * (namely, must be a multiple of 16) in the SIMD version, since it
 * refers to the address of a 128-bit integer.  In the standard C
 * version, the pointer is arbitrary.
 *
 * @param size the number of 64-bit pseudorandom integers to be
 * generated.  size must be a multiple of 2, and greater than or equal
 * to (MEXP / 128 + 1) * 2
 *
 * @note \b memalign or \b posix_memalign is available to get aligned
 * memory. Mac OSX doesn't have these functions, but \b malloc of OSX
 * returns the pointer to the aligned memory block.
 */
void fill_array64(uint64_t *array, int size) {
    assert(initialized);
    assert(idx == N32);
    assert(size % 2 == 0);
    assert(size >= N64);

    gen_rand_array((w128_t *)array, size / 2);
    idx = N32;

#if defined(BIG_ENDIAN64) && !defined(ONLY64)
    swap((w128_t *)array, size /2);
#endif
}

/**
 * This function initializes the internal state array with a 32-bit
 * integer seed.
 *
 * @param seed a 32-bit integer used as the seed.
 */
void init_gen_rand(uint32_t seed) {
    int i;

    psfmt32[idxof(0)] = seed;
    for (i = 1; i < N32; i++) {
	psfmt32[idxof(i)] = 1812433253UL * (psfmt32[idxof(i - 1)]
					    ^ (psfmt32[idxof(i - 1)] >> 30))
	    + i;
    }
    idx = N32;
    period_certification();
    initialized = 1;
}

/**
 * This function initializes the internal state array,
 * with an array of 32-bit integers used as the seeds
 * @param init_key the array of 32-bit integers, used as a seed.
 * @param key_length the length of init_key.
 */
void init_by_array(uint32_t *init_key, int key_length) {
    int i, j, count;
    uint32_t r;
    int lag;
    int mid;
    int size = N * 4;

    if (size >= 623) {
	lag = 11;
    } else if (size >= 68) {
	lag = 7;
    } else if (size >= 39) {
	lag = 5;
    } else {
	lag = 3;
    }
    mid = (size - lag) / 2;

    memset(sfmt, 0x8b, sizeof(sfmt));
    if (key_length + 1 > N32) {
	count = key_length + 1;
    } else {
	count = N32;
    }
    r = func1(psfmt32[idxof(0)] ^ psfmt32[idxof(mid)]
	      ^ psfmt32[idxof(N32 - 1)]);
    psfmt32[idxof(mid)] += r;
    r += key_length;
    psfmt32[idxof(mid + lag)] += r;
    psfmt32[idxof(0)] = r;

    count--;
    for (i = 1, j = 0; (j < count) && (j < key_length); j++) {
	r = func1(psfmt32[idxof(i)] ^ psfmt32[idxof((i + mid) % N32)]
		  ^ psfmt32[idxof((i + N32 - 1) % N32)]);
	psfmt32[idxof((i + mid) % N32)] += r;
	r += init_key[j] + i;
	psfmt32[idxof((i + mid + lag) % N32)] += r;
	psfmt32[idxof(i)] = r;
	i = (i + 1) % N32;
    }
    for (; j < count; j++) {
	r = func1(psfmt32[idxof(i)] ^ psfmt32[idxof((i + mid) % N32)]
		  ^ psfmt32[idxof((i + N32 - 1) % N32)]);
	psfmt32[idxof((i + mid) % N32)] += r;
	r += i;
	psfmt32[idxof((i + mid + lag) % N32)] += r;
	psfmt32[idxof(i)] = r;
	i = (i + 1) % N32;
    }
    for (j = 0; j < N32; j++) {
	r = func2(psfmt32[idxof(i)] + psfmt32[idxof((i + mid) % N32)]
		  + psfmt32[idxof((i + N32 - 1) % N32)]);
	psfmt32[idxof((i + mid) % N32)] ^= r;
	r -= i;
	psfmt32[idxof((i + mid + lag) % N32)] ^= r;
	psfmt32[idxof(i)] = r;
	i = (i + 1) % N32;
    }

    idx = N32;
    period_certification();
    initialized = 1;
}

/**
 *
 * This function certificate the period of 2^{MEXP}
 */
static void period_certification(void) {
    int inner = 0;
    int i, j;
    uint32_t work;

    for (i = 0; i < 4; i++)
	inner ^= psfmt32[idxof(i)] & parity[i];
    for (i = 16; i > 0; i >>= 1)
	inner ^= inner >> i;
    inner &= 1;
    /* check OK */
    if (inner == 1) {
	return;
    }
    /* check NG, and modification */
    for (i = 0; i < 4; i++) {
	work = 1;
	for (j = 0; j < 32; j++) {
	    if ((work & parity[i]) != 0) {
		psfmt32[idxof(i)] ^= work;
		return;
	    }
	    work = work << 1;
	}
    }
}

// disable little trick
#undef static

}; // end of class declaration;

} // end of namespace

#endif // SFMT_HPP
