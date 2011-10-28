/*
 * univariate_likelihodd.hpp
 *
 * The implementation of the Univariate Likelihood score L3
 *
 *  Created on: Jun 24, 2010
 *      Author: kforner
 */

#ifndef UNIVARIATE_LIKELIHOOD_HPP_
#define UNIVARIATE_LIKELIHOOD_HPP_

#include "univariate_table.hpp"
#include "matrix.hpp"

class GenotypeErrorModel;

class UnivariateLikelihood  {
public: // ========== LIFECYCLE ==========
	UnivariateLikelihood(int nb_phenotypes, double regCoeff);
	virtual ~UnivariateLikelihood();

public: // =========== ACCESSORS ==========
	int nb_phenotypes() const;
	double regCoeff() const;

public:
	// =========== METHODS ==========
	double logLikelihood(const UnivariateTable& table,
			Matrix<double> const& error_model, bool only_h1 = false);

	double logLikelihood(const vector<vector<UBYTE> >& genos_by_snp,
			const vector<UBYTE>& phenotypes,
			GenotypeErrorModel const& error_models, bool only_h1 = false);

	double logLikelihood(const vector<UnivariateTable>& tables,
			GenotypeErrorModel const& error_models, bool only_h1 = false);

	friend std::ostream& operator<<(std::ostream& o,
			UnivariateLikelihood const& u);

public: // ========== CLASS METHODS ==========
	/*
	 * Compute the frequencies of the "true" genotypes from the
	 * counts of observed ones per phenotype
	 *
	 * @param 	table		the contingency table
	 * @param   regCoeff 	regularization rate: will be used to compute the constant to add
	 * 						in each cell of the contingency table to regularize estimators. This constant will
	 * 						be = rate*total/nb_of_cells
	 *
	 * @param	[out]probs	the probs/frequencies
	 */
	static void computeUnivariateGenotypeFrequencies(const UnivariateTable& table, double regCoeff, Matrix<double>& probs, bool only_h1=false);

	/*
	 * Compute the table of likelihood of observed genotype relative to the phenotypes
	 *
	 *  The likelihood of of an observed genotype O for a phenotype s is L(O\s)= sum_over all G (P(O\G\s)*P(G\s) )
	 *
	 * @param 	error_model				the error model: P(O\G), probabilities of observed genotypes knowing true genotypes
	 * 									See fillGenotypeLikelihoodMatrixForAffymetrixModel() and fillGenotypeLikelihoodMatrixForSimpleErrorModel()
	 * 									size: [NB_GENOTYPES_WITH_ZZ][NB_GENOTYPES_NO_ZZ]
	 * @param   genotype_probabilities 	the (estimated) probabilities of true genotypes relative to the phenotype
	 *
	 * @param	[out]geno_llhs			the likelihoods per genotype and per phenotype
	 */
	static void computeUnivariateLikelihood_by_GenoPhenotype(Matrix<double> const& error_model,
			const Matrix<double>& probs, Matrix<double>& geno_llhs, bool only_h1=false);

	/*

	 * Compute the log-likelihood of observed genotype relative to the phenotypes for
	 * all patients on ONE SNP
	 *
	 *  The likelihood of of an observed genotype O for a phenotype s is L(O\s)= sum_over all G (P(O\G\s)*P(G\s) )
	 *  So The likelihood for all samples is the product over all sample genotypes O_k product(L(O_k\s)
	 *
	 * @param   contingency_table		the table as output by fillGenotypeContingencyTableOneSNPBinaryLabel
	 * @param 	genotype_likelihoods	the likelihoods per genotype and per phenotype (see computeUnivariateLikelihood_by_GenoPhenotype() )
	 *

	 * @return  the log_likelihood
	 */

	static double computeUnivariateLogLikelihoodForAllSamples_H1(
			const UnivariateTable& table,
			const Matrix<double>& geno_llhs);

	/*
	 * Function: computeUnivariateLogLikelihoodForAllSamples_H1
	 *
	 * Compute the log-likelihood of observed genotype without considering the phenotype
	 *
	 * @param   contingency_table           the table as output by fillGenotypeContingencyTableOneSNPBinaryLabel
	 * @param       genotype_likelihoods    the likelihoods per genotype and per phenotype (see computeUnivariateLikelihood_by_GenoPhenotype() )
	 *

	 * @return  the log_likelihood
	 */

	static double computeUnivariateLogLikelihoodForAllSamples_H0(
			const UnivariateTable& table,
			const Matrix<double>& geno_llhs);


	static double computeUnivariateLogLikelihoodForAllSamples(
			const UnivariateTable& table,
			const Matrix<double>& geno_llhs, bool onlyH1 = false);

	/*
	 * Compute the log-likelihood ratio for all samples
	 *
	 *
	 * @param   contingency_table           the table as output by fillGenotypeContingencyTableOneSNPBinaryLabel
	 * @param       genotype_likelihoods    the likelihoods per genotype and per phenotype (see computeUnivariateLikelihood_by_GenoPhenotype() )
	 *

	 * @return  the log_likelihood
	 */

	static double computeUnivariateLogLikelihoodRatioForAllSamples(
			const UnivariateTable& table,
			const Matrix<double>& geno_llhs);


private: // ========== INSTANCE DATA ==========
	double _reg_coeff;
	UnivariateTable _table;
	Matrix<double> _probs;
	Matrix<double> _geno_llhs;
};

// ========== UnivariateLikelihood INLINE IMPLEMENTATION  ==========
inline int UnivariateLikelihood::nb_phenotypes() const {
	return _table.nb_phenotypes();
}

inline double UnivariateLikelihood::regCoeff() const {
	return _reg_coeff;
}


inline double UnivariateLikelihood::computeUnivariateLogLikelihoodForAllSamples(
		const UnivariateTable & table, const Matrix<double> & geno_llhs,
		bool onlyH1) {
	return onlyH1 ? computeUnivariateLogLikelihoodForAllSamples_H1(table,
			geno_llhs)
			: computeUnivariateLogLikelihoodRatioForAllSamples(table,
					geno_llhs);
}


inline double UnivariateLikelihood::computeUnivariateLogLikelihoodRatioForAllSamples(
		const UnivariateTable & table, const Matrix<double> & geno_llhs) {
	double h1 =
			computeUnivariateLogLikelihoodForAllSamples_H1(table, geno_llhs);
	double h0 =
			computeUnivariateLogLikelihoodForAllSamples_H0(table, geno_llhs);
//	cerr << h1 << " " << h0 << endl;
	return h1 - h0;
}


#endif /* UNIVARIATE_LIKELIHOOD_HPP_ */
