
/*
 * sum_scores.cpp
 *
 *  Created on: Jul 7, 2010
 *      Author: kforner
 */

#include "sum_scores.hpp"
#include "misc.hpp"


void AllelicSumScore::precompute_allelic_expected_frequencies(
		const UnivariateTable& table, double** expected)
{
	const int nb_true_phenotypes = table.nb_phenotypes();
	const int h0_phenotype = table.h0_phenotype();

	int total_a = 2*table[h0_phenotype][HOMOZYGOTE1] + table[h0_phenotype][HETEROZYGOTE];
	int total_A = 2*table[h0_phenotype][HOMOZYGOTE2] + table[h0_phenotype][HETEROZYGOTE];
	double total = total_a + total_A;

	if ( total == 0 ) {
		for (int k = 0; k < nb_true_phenotypes; ++k)
			expected[k][ALLELE1] = expected[k][ALLELE2] = 0;
	} else {

		for (int k = 0; k < nb_true_phenotypes; ++k) {
			int row_sum = 0;
			for (int i = 0; i < NB_GENOTYPES_NO_ZZ; ++i)
				row_sum += table[k][i];
			// needs to double to get the count of alleles
			row_sum += row_sum;

			expected[k][ALLELE1] = (total_a*row_sum)/total;
			expected[k][ALLELE2] = (total_A*row_sum)/total;
		}
	}

}

double AllelicSumScore::allelic_pearson_score(const UnivariateTable& table,
		const double*const* expected)
{

	double score = 0;
	const int nb_true_phenotypes = table.nb_phenotypes();
	for (int k = 0; k < nb_true_phenotypes; ++k) {
		const int* row = table[k];
		double e1 = expected[k][ALLELE1];
		double e2 = expected[k][ALLELE2];
		int nb_het = row[HETEROZYGOTE];

		if ( e1 > 0 ) {
			int nb_a = 2*row[HOMOZYGOTE1] + nb_het;
			double s = nb_a - e1;
			score += (s*s)/e1;
		}

		if ( e2 > 0 ) {
			int nb_A = 2*row[HOMOZYGOTE2] + nb_het;
			double s = nb_A - e2;
			score += (s*s)/e2;
		}
	}

	return score;
}


void GenotypicSumScore::precompute_genotypic_expected_frequencies(
		const UnivariateTable& table, double** expected)
{

	double total = table.total(false);
	const int nb_true_phenotypes = table.nb_phenotypes();
	const int h0_phenotype = table.h0_phenotype();

	if ( total == 0 ) {
		for (int k = 0; k < nb_true_phenotypes; ++k)
			for (int i = 0; i < NB_GENOTYPES_NO_ZZ; ++i)
				expected[k][i] = 0;
	} else {
		for (int k = 0; k < nb_true_phenotypes; ++k) {
			int row_sum = 0;
			for (int i = 0; i < NB_GENOTYPES_NO_ZZ; ++i)
				row_sum += table(k, i);
			for (int i = 0; i < NB_GENOTYPES_NO_ZZ; ++i)
				expected[k][i] = (table(h0_phenotype, i) * row_sum) / total;
		}
	}
}

double GenotypicSumScore::genotypic_pearson_score(const UnivariateTable& table, const double*const* expected)
{

	double score = 0;
	const int nb_true_phenotypes = table.nb_phenotypes();
	for (int k = 0; k < nb_true_phenotypes; ++k) {
		for (int i = 0; i < NB_GENOTYPES_NO_ZZ; ++i) {
			double e = expected[k][i];
			if ( e > 0 ) {
				double s = table[k][i] - e;
				score += (s*s)/e;
			}
		}
	}

	return score;
}

