/*
 * sum_scores.hpp
 *
 * Define the AllelicSumScore and GenotypicSumScore statistics.
 *
 *  Created on: Jul 7, 2010
 *      Author: kforner
 */

#ifndef SUM_SCORES_HPP_
#define SUM_SCORES_HPP_

//#include "pvalues.hpp"
#include "univariate_table.hpp"
#include "matrix3d.hpp"

class AllelicSumScore  {

public: // ========== CLASS METHODS ==========
	static void precompute_allelic_expected_frequencies(const UnivariateTable& table, double** expected);
	static double allelic_pearson_score(const UnivariateTable& table, const double*const* expected);

public: // ========== CONSTANTS ==========
	static const int NB_ALLELES_NO_Z = 2;
	static const int ALLELE1 = 0;
	static const int ALLELE2 = 1;
};

class GenotypicSumScore  {
public: // ========== CLASS METHODS ==========
	static void precompute_genotypic_expected_frequencies(const UnivariateTable& table, double** expected);
	static double genotypic_pearson_score(const UnivariateTable& table, const double*const* expected);
};

#endif /* SUM_SCORES_HPP_ */
