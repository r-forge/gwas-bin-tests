/*
 * types.hpp
 *
 *  Created on: Jun 14, 2010
 *      Author: kforner
 */

#ifndef TYPES_HPP_
#define TYPES_HPP_

#include <cfloat>
#include <stdint.h>

typedef int64_t INT64;

/**Integer type for data set.
*/
typedef unsigned char UBYTE;

/**Floating point number
*/
typedef double UDOUBLE;
/** Max value for floating point number
*/


typedef unsigned long ULONG;

#define UDBL_MAX ((UDOUBLE)DBL_MAX)
 /*const UDOUBLE UDBL_MAX = DBL_MAX; */

#define HOMOZYGOTE1 	0
#define HETEROZYGOTE 	1
#define HOMOZYGOTE2 	2
#define MISSING			3

#define NB_GENOTYPES_NO_ZZ	3
#define NB_GENOTYPES_WITH_ZZ 4
/*#define NB_PHENOTYPES 2*/

#define DEBUG_MODE 1

#endif /* TYPES_HPP_ */
