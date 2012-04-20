//#define MEXP 607

//namespace NAMESPACE {

// brings all static global variables and inline static funtions


//#define static
#include "SFMT.c"
//#undef static

	class RNG {
private: // ========== DECLARATION OF INSTANCE STATE VARIABLES =========

		// ===== instance variable ====
		// These state variables come from the "FILE GLOBAL VARIABLES" section
		// in SFMT.c: take all variables without the "static" keyword
		// and put the definitions in the init() method
		// because it is not allowed  in C++ to define vars inside a class

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



	public: // ========== LIFECYCLE
		/**
		 * init the rng using SFMT C function gen_rand_all() :
		 * This function fills the internal state array with pseudorandom
		 * integers.
		 */
		RNG() { init(); gen_rand_all(); }

		/**
		 * init the rng with a seed
		 * @see set_seed
		 */
		RNG(uint32_t seed) { init(); set_seed(seed); }

		/**
		 * init the rng with an array of seeds
		 * @see set_seeds
		 */
		RNG(uint32_t *init_key, int key_length) {
			init();
			set_seeds(init_key, key_length);
		}



	public: // ===== PUBLIC INTERFACE =====
		/**
		 * This function generates and returns 32-bit pseudorandom number.
		 *
		 */
		uint32_t rand32() {
			// compute state sum
			uint32_t sum1 = 0;
		    for (int i = 0; i < N32; i++) {
		    	sum1 += psfmt32[idxof(i)];
		    }

			uint32_t r = gen_rand32();

			uint32_t sum2 = 0;
		    for (int i = 0; i < N32; i++) {
		    	sum2 += psfmt32[idxof(i)];
		    }
			fprintf(stderr, "idx=%i N32=%i sum1=%i sum2=%i r=%i\n", idx, N32, sum1, sum2, r);
			return r;
		}


		/**
		 * set the seed. Can be used at any time to reset the seed
		 *
		 * Call the init_gen_rand() SFMT function:
		 *
		 * This function initializes the internal state array with a 32-bit
		 * integer seed.
		 * @param seed a 32-bit integer used as the seed.
		 */
		void set_seed(uint32_t seed) { init_gen_rand(seed); }

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

		std::string id() const {
			return get_idstring();
		}

	protected:
		/**
		 * Init state variables: initialization as copy/paste
		 * of section "FILE GLOBAL VARIABLES" in SFMT.c
		 *
		 */
		 void init() {
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
		}

	private: // ==== original inline static C SFMT functions =====
		// these are the function that use state global variables
		// because we have instance variables with the same name
		// we can use the exact same code
		// so there functions are just copy/pasted here

		// Functions: gen_rand32, gen_rand64, fill_array32, fill_array64
		// init_gen_rand, init_by_array, gen_rand_all, gen_rand_array

// little trick: remove the static statements
#define static

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



//};

