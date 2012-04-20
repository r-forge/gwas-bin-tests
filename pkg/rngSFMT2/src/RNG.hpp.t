// not meant to be included directly
#include <inttypes.h>
#include <string>
class RNG {
	public: // ========== LIFECYCLE
		/**
		 * init the rng using SFMT C function gen_rand_all() :
		 * This function fills the internal state array with pseudorandom
		 * integers.
		 */
		RNG();

		/**
		 * init the rng with a seed
		 * @see set_seed
		 */
		RNG(uint32_t seed);

		/**
		 * init the rng with an array of seeds
		 * @see set_seeds
		 */
		RNG(uint32_t *init_key, int key_length);

	public: // ===== PUBLIC INTERFACE =====
		/**
		 * This function generates and returns 32-bit pseudorandom number.
		 *
		 */
		uint32_t rand32();

		/**
		 * set the seed. Can be used at any time to reset the seed
		 *
		 * Call the init_gen_rand() SFMT function:
		 *
		 * This function initializes the internal state array with a 32-bit
		 * integer seed.
		 * @param seed a 32-bit integer used as the seed.
		 */
		void set_seed(uint32_t seed);

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
		void set_seeds(uint32_t *init_key, int key_length);

		std::string id() const;
};
