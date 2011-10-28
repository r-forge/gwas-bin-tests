/*
 * univariate_likelihood.cpp
 *
 *  Created on: Jun 24, 2010
 *      Author: kforner
 */

#include "univariate_likelihood.hpp"
#include <cmath>
#include "genotype_error_model.hpp"

UnivariateLikelihood::UnivariateLikelihood(int nb_phenotypes, double regCoeff)
 : _reg_coeff(regCoeff), _table(nb_phenotypes),
   _probs(nb_phenotypes+1, NB_GENOTYPES_NO_ZZ),
   _geno_llhs(nb_phenotypes+1, NB_GENOTYPES_WITH_ZZ)
{
}

UnivariateLikelihood::~UnivariateLikelihood() {}


std::ostream& operator<<(std::ostream& o,UnivariateLikelihood const& u) {
	o << "UnivariateLikelihood:\n";
	o << "	regCoeff=" << u.regCoeff() << "\n";
	o << "table:" << u._table;
	o << "probs:" << u._probs;
	o << "likelihoods by genotype:" << u._geno_llhs;
	return o;
}

double UnivariateLikelihood::logLikelihood(const UnivariateTable& table,
		Matrix<double> const& error_model, bool only_h1) {
	// for H1 only, no need to compute the table for H0

	computeUnivariateGenotypeFrequencies(table,  regCoeff(), _probs, only_h1);
	computeUnivariateLikelihood_by_GenoPhenotype(error_model, _probs, _geno_llhs, only_h1);

	return computeUnivariateLogLikelihoodForAllSamples(table, _geno_llhs, only_h1);
}

double UnivariateLikelihood::logLikelihood(
		const vector<vector<UBYTE> >& genos_by_snp,
		const vector<UBYTE>& phenotypes,
		GenotypeErrorModel const& error_models, bool only_h1) {

	double dataset_likelihood = 0;
	const int nb_snps = SIZE(genos_by_snp);
	for (int snp = 0; snp < nb_snps; ++snp) {
		_table.fill( genos_by_snp[snp], phenotypes);
		dataset_likelihood += logLikelihood(_table, error_models.getErrorModelForSnp(snp), only_h1);
	}
	return dataset_likelihood;
}

double UnivariateLikelihood::logLikelihood(const vector<UnivariateTable>& tables,
		GenotypeErrorModel const& error_models, bool only_h1) {

	double dataset_likelihood = 0;
	const int nb_tables = SIZE(tables);
	for (int i = 0; i < nb_tables; ++i)
		dataset_likelihood += logLikelihood(tables[i], error_models.getErrorModelForSnp(i), only_h1);

	return dataset_likelihood;
}


void UnivariateLikelihood::computeUnivariateGenotypeFrequencies(
		const UnivariateTable & table, double regCoeff, Matrix<double> & probs, bool only_h1) {

	assert( table.nb_phenotypes() + 1 <= probs.nb_rows() );
	assert( probs.nb_cols() >= NB_GENOTYPES_NO_ZZ );

	const int nb = only_h1 ? table.nb_phenotypes() : table.nb_phenotypes_with_h0();

	// the regularization constant for one phenotype is: expected nb by cell times the regCoeff
	const double reg_table = regCoeff / (double)NB_GENOTYPES_NO_ZZ;

	const int* row_table = 0;
	double* row_probs;
	double sum = 0;

	for (int j = 0; j < nb; ++j) {
		sum = 0;
		row_table = table.getRow(j);
		row_probs = probs[j];
		for (int i = 0; i < NB_GENOTYPES_NO_ZZ; ++i)
			sum += row_table[i];

		double reg_constant = sum * reg_table;
		sum += NB_GENOTYPES_NO_ZZ*reg_constant;
		for (int i = 0; i < NB_GENOTYPES_NO_ZZ; ++i)
			row_probs[i] = (row_table[i]+reg_constant) / sum;
	}
}


void UnivariateLikelihood::computeUnivariateLikelihood_by_GenoPhenotype(
		Matrix<double> const& error_model,
		const Matrix<double>& probs,
		Matrix<double>& geno_llhs, bool only_h1) {

	double sum;
	const double* row_probs;
	double* row_results;
	const int nb = only_h1 ? probs.nb_rows() - 1 : probs.nb_rows();
	for (int S = 0; S < nb; ++S)
	{
		row_probs = probs[S];
		row_results = geno_llhs[S];
		for (int O = 0; O < NB_GENOTYPES_WITH_ZZ; ++O)
		{
			sum = 0;
			for (int G = 0; G < NB_GENOTYPES_NO_ZZ; ++G)
				sum += error_model[O][G] * row_probs[G];
			row_results[O] = sum;

		}
	}
}


double UnivariateLikelihood::computeUnivariateLogLikelihoodForAllSamples_H1(
		const UnivariateTable& table,
		const Matrix<double>& geno_llhs) {

	const int* row_table;
	const double* row_gl;
	double llh = 0;
	int nb = 0;
	double lh_k_g = 0;
	const int K = table.nb_phenotypes();
	for (int k = 0; k < K; ++k) {
		row_table = table.getRow(k);
		row_gl = geno_llhs[k];
		for (int g = 0; g < NB_GENOTYPES_WITH_ZZ; ++g) {
			nb = row_table[g];
			if (nb > 0) {
				lh_k_g = row_gl[g];
				if ( lh_k_g > 0 )
					llh += nb * log(lh_k_g);
			}
		}
	}
	return llh;
}

double UnivariateLikelihood::computeUnivariateLogLikelihoodForAllSamples_H0(
		const UnivariateTable& table,
		const Matrix<double>& geno_llhs) {

	const int* row_table = table.getH0Row();
	const double* row_gl = geno_llhs[table.h0_phenotype()];

	int nb = 0;
	double llh = 0;
	double lh_g = 0;
	for (int g = 0; g < NB_GENOTYPES_WITH_ZZ; ++g) {
		nb = row_table[g];
		if (nb > 0) {
			lh_g = row_gl[g];
			if ( lh_g > 0 )
				llh += nb * log(lh_g);
		}
	}

	return llh;
}

