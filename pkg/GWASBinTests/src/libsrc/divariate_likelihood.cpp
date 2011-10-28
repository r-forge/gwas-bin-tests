/*
 * divariate_likelihood.cpp
 *
 *  Created on: Jul 6, 2010
 *      Author: kforner
 */

#include <cmath>
#include "divariate_likelihood.hpp"

DivariateLikelihood::DivariateLikelihood(int nb_phenotypes, double regCoeff)
	: _reg_coeff(regCoeff),
	  _probs(nb_phenotypes+1, NB_GENOTYPES_NO_ZZ, NB_GENOTYPES_NO_ZZ),
	  _geno_llhs(nb_phenotypes+1, NB_GENOTYPES_WITH_ZZ, NB_GENOTYPES_WITH_ZZ)
{}

DivariateLikelihood::~DivariateLikelihood() {}


// =========== METHODS ==========
double DivariateLikelihood::logLikelihood(const DivariateTable& table, Matrix<
		double> const& error_model1, Matrix<double> const& error_model2,
		 bool only_h1)
{
	computeDivariateGenotypeFrequencies(table, regCoeff(), _probs, only_h1);
	computeDivariateLikelihood_by_GenosPhenotype(error_model1, error_model2, _probs,
			_geno_llhs, only_h1);
	return computeDivariateLogLikelihoodForAllSamples(table, _geno_llhs,
			only_h1);
}

double DivariateLikelihood::computeDivariateLogLikelihoodForAllSamples(
		const DivariateTable& table, const Matrix3D<double>& geno_llhs,
		bool onlyH1)
{
	return onlyH1 ?
			computeDivariateLogLikelihoodForAllSamples_H1(table, geno_llhs)
		: computeDivariateLogLikelihoodForAllSamples_H0(table, geno_llhs);
}

double DivariateLikelihood::computeDivariateLogLikelihoodForAllSamples_H1(
		const DivariateTable& table, const Matrix3D<double>& geno_llhs)
{
	double llh = 0;
	double lh_k_o1_o2 = 0;
	int nb = 0;
	const int nb_true_phenotypes = table.nb_phenotypes();
	for (int k = 0; k < nb_true_phenotypes; ++k)
	{
		const int*const* row_table = table.getMatrix(k);
		const double*const* row_gl = geno_llhs[k];
		for (int o1 = 0; o1 < NB_GENOTYPES_WITH_ZZ; ++o1)
		{
			const int* col_table = row_table[o1];
			const double* col_gl = row_gl[o1];
			for (int o2 = 0; o2 < NB_GENOTYPES_WITH_ZZ; ++o2)
			{
				nb = col_table[o2];

				if (nb > 0) {
					lh_k_o1_o2 = col_gl[o2];
					if (lh_k_o1_o2 > 0)
					llh += nb * log(lh_k_o1_o2);
				}

			}
		}
	}

	return llh;
}

double DivariateLikelihood::computeDivariateLogLikelihoodForAllSamples_H0(
		const DivariateTable& table, const Matrix3D<double>& geno_llhs)
{

	double llh_h0 = 0;
	double lh_o1_o2 = 0;
	int nb = 0;

	const double*const* row_gl = geno_llhs[table.h0_phenotype()];
	const int*const* row_table = table.getMatrix( table.h0_phenotype() );
	for (int o1 = 0; o1 < NB_GENOTYPES_WITH_ZZ; ++o1) {
		const int* col_table = row_table[o1];
		const double* col_gl = row_gl[o1];
		for (int o2 = 0; o2 < NB_GENOTYPES_WITH_ZZ; ++o2) {
			nb = col_table[o2];
			if (nb > 0) {
				lh_o1_o2 = col_gl[o2];
				if (lh_o1_o2 > 0)
					llh_h0 += nb * log(lh_o1_o2);
			}
		}
	}

	return llh_h0;
}

void DivariateLikelihood::computeDivariateGenotypeFrequencies(
		const DivariateTable& table, double regCoeff, Matrix3D<double>& probs,
		bool only_h1)
{
	assert( table.nb_phenotypes() + 1 <= probs.nx() );
	assert( probs.ny() >= NB_GENOTYPES_NO_ZZ );
	assert( probs.nz() >= NB_GENOTYPES_NO_ZZ );

	const int NB_CELLS = NB_GENOTYPES_NO_ZZ*NB_GENOTYPES_NO_ZZ;
	const double reg_table = regCoeff / (double)NB_CELLS;
	const int nb = only_h1 ? table.nb_phenotypes() : table.nb_phenotypes_with_h0();

	int const*const* t_k = 0;
	double** p_k = 0;
	double* p_k_i = 0;
	int const* t_k_i = 0;
	double sum = 0;

	for (int k = 0; k < nb; ++k)
	{
		t_k = table.getMatrix(k);
		p_k = probs[k];

		// compute the sum of the table for phenotype k
		sum = 0;
		for (int i = 0; i < NB_GENOTYPES_NO_ZZ; ++i) {
			t_k_i = t_k[i];
			for (int j = 0; j < NB_GENOTYPES_NO_ZZ; ++j)
				sum += t_k_i[j];
		}

		//then compute the regularization constant
		double reg_constant = sum * reg_table;
		sum += NB_CELLS*reg_constant;

		// now we can compute the regularized estimates

		for (int i = 0; i < NB_GENOTYPES_NO_ZZ; ++i) {
			t_k_i = t_k[i];
			p_k_i = p_k[i];
			for (int j = 0; j < NB_GENOTYPES_NO_ZZ; ++j)
				p_k_i[j] = (t_k_i[j] + reg_constant) / sum;
		}
	}
}

void DivariateLikelihood::computeDivariateLikelihood_by_GenosPhenotype(Matrix<
		double> const& error_model1, Matrix<double> const& error_model2,
		const Matrix3D<double>& probs, Matrix3D<double>& results,
		bool only_h1)
{
	double sum;
	double P_O1_G1;

	const int nb = only_h1 ? probs.nx() - 1 : probs.nx();
	for (int k = 0; k < nb; ++k) {
		const double*const* genotype_div_probs = probs[k];
		double** L = results[k];
		for (int O1 = 0; O1 < NB_GENOTYPES_WITH_ZZ; ++O1)
		{
			double* L_O1 = L[O1];
			const double* P_O1 = error_model1[O1]; /* P[O1\.] */
			for (int O2 = 0; O2 < NB_GENOTYPES_WITH_ZZ; ++O2)
			{
				const double* P_O2 = error_model2[O2]; /* P[O2\.] */
				sum = 0;
				for (int G1 = 0; G1 < NB_GENOTYPES_NO_ZZ; ++G1)
				{
					P_O1_G1 = P_O1[G1]; /* P[O1\G1] */
					const double* P_G1 = genotype_div_probs[G1];
					for (int G2 = 0; G2 < NB_GENOTYPES_NO_ZZ; ++G2)
						sum += P_O2[G2]*P_G1[G2] * P_O1_G1; /* P[O2\G2]*P[G1,G2]*P[O1\G1] */
				}
				L_O1[O2] = sum; /* L(O1,O2) */
			}
		}

	}
}


