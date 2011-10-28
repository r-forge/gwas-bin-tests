/*
 * divariate_likelihood.hpp
 *
 *  Created on: Jul 6, 2010
 *      Author: kforner
 */

#ifndef DIVARIATE_LIKELIHOOD_HPP_
#define DIVARIATE_LIKELIHOOD_HPP_

#include "matrix3d.hpp"
#include "matrix.hpp"

#include "genotype_error_model.hpp"
#include "divariate_table.hpp"

class DivariateLikelihood {
public: // ========== LIFECYCLE ==========
	DivariateLikelihood(int nb_phenotypes, double regCoeff);
	virtual ~DivariateLikelihood();

public: // =========== ACCESSORS ==========
	int nb_phenotypes() const;
	double regCoeff() const;

public: // =========== METHODS ==========
	double logLikelihood(const DivariateTable& table,
			Matrix<double> const& error_model1, Matrix<double> const& error_model2, bool only_h1 = false);


public: // ========== CLASS METHODS ==========
	/*

	 * Compute the estimations of the joint  probabilities P(G_j-1, Gj,s) of the "true" genotypes
	 * based on the counts of the observed genotypes (Somewhat related to formula (15) in the article)
	 *
	 * The ouput matrix is probs[Genotype for SNP#1][Genotype for SNP#2][Phenotype]
	 *  i.e P[G2,G1\S]
	 *
	 * Caution:
	 *
	 *  the first dimension is for genotypes of SNP #1 !!!
	 *  the second dimension n is for genotypes of SNP #2
	 *  	and the 3rd for the phenotype
	 *
	 * @param	genotypes1	the genotypes for the SNP #1
	 * @param	genotypes2	the genotypes for the SNP #2
	 * @param	phenotype	the binary phenotype
	 * @param   size 		the size of the genotypes and phenotypes array
	 *
	 *
	 * @param 	table		the contingency table
	 * @param   regCoeff 	regularization constant added in each cell of the contingency table to regularize estimators.
	 *
	 * @param	[out]probs	the probs/frequencies
	 */
	static void computeDivariateGenotypeFrequencies(const DivariateTable& table,  double regCoeff, Matrix3D<double>& probs, bool only_h1=false);

	/*
	 * Compute the table of likelihood of a pair of observed genotypes for a SNP 1 and its successor 2
	 *
	 * The likelihood of an observed pair of genotypes (O1, O2) for a phenotype s is:
	 *  L(O1,O2\s) = sum over all (G1,G2) of P[O2\G2]*P[G2,G1]*P[O1\G1]
	 *
	 * @param 	error_model1			the error model for SNP1: P(O1\G1), probabilities of observed genotypes knowing true genotypes
	 * 									See fillGenotypeLikelihoodMatrixForAffymetrixModel() and fillGenotypeLikelihoodMatrixForSimpleErrorModel()
	 * 									size: [NB_GENOTYPES_WITH_ZZ][NB_GENOTYPES_NO_ZZ]
	 *
	 * @param 	error_model2			the error model for SNP2
	 *
	 * @param   genotype_div_probs		the (estimated) joint probabilities of P[G2,G1\S]
	 *
	 * @param	[out]results			the likelihoods per pair of genotypes genotype and per phenotype R[O1][O2][S]
	 */
	static void computeDivariateLikelihood_by_GenosPhenotype(Matrix<double> const& error_model1, Matrix<double> const& error_model2, const Matrix3D<double>& probs, Matrix3D<double>& geno_llhs, bool only_h1=false);

	/*
	 * compute the ration LikelihoodForAllSamples_H1/LikelihoodForAllSamples_H1
	 *
	 * @param   contingency_table		the table as output by fillGenotypeContingencyTableTwoSNPsBinaryLabel
	 * @param 	genotypes_likelihoods	the likelihoods per pair of genotype and per phenotype (see computeGenotypesDivariateLikelihood )
	 *

	 * @return  the log likelihood ratio
	 */
	static double computeDivariateLogLikelihoodRatioForAllSamples(const DivariateTable& table, const Matrix3D<double>& geno_llhs);

	static double computeDivariateLogLikelihoodForAllSamples(const DivariateTable& table, const Matrix3D<double>& geno_llhs, bool onlyH1 = false);

	/*
	 * Compute the log-likelihood of pairs of observed genotypes relative to the phenotypes for
	 * all patients on a pair of SNPs
	 *
	 *  The likelihood of an observed pair of genotypes O1,O2 for a phenotype s is L(O1,O2\s)
	 *  So The likelihood for all samples is the product over all samples pairs of genotypes (O1k,O2k)
	 *  of L(O1k,O2k\S)
	 *
	 * @param   contingency_table		the table as output by fillGenotypeContingencyTableTwoSNPsBinaryLabel
	 * @param 	genotypes_likelihoods	the likelihoods per pair of genotype and per phenotype (see computeGenotypesDivariateLikelihood )
	 *

	 * @return  the log_likelihood
	 */

	static double computeDivariateLogLikelihoodForAllSamples_H1(const DivariateTable& table, const Matrix3D<double>& geno_llhs);

	/*
	 * Compute the log-likelihood of pairs of observed genotypes under H0 i.e. without
	 * taking into account the phenotype
	 *
	 * @param   contingency_table		the table as output by fillGenotypeContingencyTableTwoSNPsBinaryLabel
	 * @param 	genotypes_likelihoods	the likelihoods per pair of genotype and per phenotype (see computeGenotypesDivariateLikelihood )
	 *

	 * @return  the log_likelihood
	 */
	static double computeDivariateLogLikelihoodForAllSamples_H0(const DivariateTable& table, const Matrix3D<double>& geno_llhs);


private: // ========== INSTANCE DATA ==========
	double _reg_coeff;
	Matrix3D<double> _probs;
	Matrix3D<double> _geno_llhs;
};

// =========== ACCESSORS ==========
inline int DivariateLikelihood::nb_phenotypes() const {
	return  _probs.nx() - 1;
}

inline double DivariateLikelihood::regCoeff() const {
	return _reg_coeff;
}


#endif /* DIVARIATE_LIKELIHOOD_HPP_ */
