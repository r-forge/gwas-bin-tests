/*
 * genotype_error_model.hpp
 *
 *  Created on: Jun 18, 2010
 *      Author: kforner
 */

#ifndef GENOTYPE_ERROR_MODEL_HPP_
#define GENOTYPE_ERROR_MODEL_HPP_


#include "matrix.hpp"

class GenotypeErrorModel {
public: // ========== LIFECYCLE ==========
	explicit GenotypeErrorModel(const vector<vector<UBYTE> >& data);
public: // =========== ACCESSORS ==========
	int getNbSnps() const;
	Matrix<double> const& getErrorModelForSnp(int snp) const;
	/*
	 * return the probability for a given snp to have genotype G
	 * if the observed genotype ie O: P(O_j \ G_j) for snp j
	 */
	double getProbability(int snp, int O, int G) const;
public: // =========== METHODS ==========
	friend std::ostream& operator<< (std::ostream& o, GenotypeErrorModel const& GenotypeErrorModel);

public: // =========== CLASS METHODS ==========

	/*
	 * Given some genotypes, count the numbers of:
	 * 		- first homozygotous (aa)
	 * 		- heterozygotous (aA)
	 * 		- second homozygotous (AA)
	 * 		- missing (ZZ)
	 *
	 * N.B: it is not exactly the same as the legacy computeSNPHomoHeteroTable()
	 *
	 *
	 * @param genotypes				genotypes
	 * @param [out]counts			the genotype counts to fill
	 *
	 *
	 */
	static void count_genotypes(vector<UBYTE> const& genotypes, int counts[NB_GENOTYPES_WITH_ZZ]);

protected:
	vector< Matrix<double> >& getMatrix();
	vector< Matrix<double> > const& getMatrix() const;

private: // ========== INSTANCE DATA ==========
	vector< Matrix<double> > _matrices;
};

class SimpleGenotypeErrorModel : public GenotypeErrorModel {
public: // ========== LIFECYCLE ==========
	explicit SimpleGenotypeErrorModel(const vector<vector<UBYTE> >& data, double max_error_rate=DEFAULT_ERROR_RATE);
public: // =========== ACCESSORS ==========
	double getMaxErrorRate() const;
//	double getBeta() const;
public: // =========== CLASS METHODS ==========
	/*
	 *  Fill the genotype likelihood matrix according to the
	 *  error Simple model with as many errors on the heterozygotous calls and on homozygotous ones
	 Corresponding matrix (equivalent of Eq. 15) :
	 (1-beta)*(1-alpha)	(1-beta)*alpha/2	(1-beta)*alpha/2
	 (1-beta)*alpha/2		(1-beta)*(1-alpha)	(1-beta)*alpha/2
	 (1-beta)*alpha/2		(1-beta)*alpha/2	(1-beta)*(1-alpha)
	 beta				beta				beta

	 *
	 * @param 	alpha			the error rate
	 * @param 	beta			the missing rate
	 * @param 					the matrix to fill, must be allocated of size [4][3] where
	 *  						n is the number of snps == (genotypes->nrows)
	 *
	 */
	static void fillGenotypeLikelihoodMatrixForSimpleErrorModel(double alpha, double beta, Matrix<double>& P);

public: // ========== CLASS CONSTANTS ==========
	static const double DEFAULT_ERROR_RATE; // = 0.05;
private: // ========== INSTANCE DATA ==========
	double _maxErrorRate;
//	double _beta;
};

class AffymetrixGenotypeErrorModel : public GenotypeErrorModel {
public: // ========== LIFECYCLE ==========
	explicit AffymetrixGenotypeErrorModel(const vector<vector<UBYTE> >& data, double max_error_rate=SimpleGenotypeErrorModel::DEFAULT_ERROR_RATE);
public: // =========== ACCESSORS ==========
	double getMaxErrorRate() const;

public: // =========== CLASS METHODS ==========
	/*
	 * When there is no solution for beta the missing rate in the affymetrix error model
	 *  (see solveAffymetrixMissingRate() ), that should happen only for SNPs with low diversity
	 *  the a-priori error rate is decreased to the maximum possible error-rate obtained for
	 *  the unrealistic value of P(O=Aa)=0
	 *
	 *  The corresponding missing rate is then P(O=missing data).
	 *
	 *  This function computes those two values
	 *
	 * @param counts			the genotype counts
	 * @param [out]				alpha
	 * @param [out]				beta
	 *
	 */
	static void fixAffymetrixRatesWhenNoSolution(int counts[NB_GENOTYPES_WITH_ZZ],
			double& alpha_out, double& beta_out,  bool debug=true);

	/*
	 *
	 *  Compute the missing rate according to the model described in p8 (eq 17) of
	 http://velblod.videolectures.net/2007/pascal/msht07_paris/omont_nicolas/article10.pdf .
	 The basic idea was that the original calling algorithm from Affymetrix tended to have a higher proportion
	 of missing values for heterozygotous. Therefore, a missing value is more likely to be an heterozygotous.
	 The other idea is that it is very unlikely to mistake an homozygotous for the other,
	 while it is possible to mistake one homozygotous for an heterozygotous.
	 *
	 * Given alpha, the error rate which :
	 * "is estimated during external comparison
	 of Affymetrix GeneChip°R human mapping
	 100K and other technologies. In this study, the error
	 rate is chosen to be a = 0:05 (cf. Affymetrix
	 GeneChip) "
	 * ,this function computes the missing rate beta by solving a second degree equation.
	 *  If this equation has no solution, see the
	 *  	fixAffymetrixRatesWhenNoSolution()
	 *
	 *
	 * N.B: it is not exactly the same as the legacy computeSNPHomoHeteroTable()
	 *
	 *
	 * @param alpha				the estimated error rate (alpha) for the dataset. Usually 5%.
	 * @param counts			the genotype counts
	 *
	 * @return the calculate missing rate (>=0) or
	 * 		-1 (<0) if no root to the equation
	 *      -2 if no suitable root (i.e no root that meets the rules fo the system
	 *
	 */
	static double solveAffymetrixMissingRate(double alpha, int counts[NB_GENOTYPES_WITH_ZZ]);
	/*
	 *
	 *  Fill the genotype likelihood matrix according to the affymetrix error model
	 *  and the error and missing rates.
	 *  This model tries to correct the bias on missing values that result in a lack of heterozygotous
	 *  The model is decribed in Eq. 15 of article10.pdf
	 *
	 *  this matrix is interpreted as follows P[i][j] is the likelihood of genotype j
	 *  for observed genotype i (using the constants HOMOZYGOTE1, HOMOZYGOTE2, HETEROZYGOTE, MISSING)
	 *  i.e P(O=i\G=j)
	 *
	 * See solveAffymetrixMissingRate() and fixAffymetrixRatesWhenNoSolution() for how
	 * to get the rates.
	 *
	 * @param 	alpha				the error rate
	 * @param 	beta				the missing rate
	 * @param	P					the matrix to fill, must be allocated of size [4][3]

	 *
	 */
	static void fillGenotypeLikelihoodMatrixForAffymetrixModel(double alpha, double beta, Matrix<double>& P);

public: // ========== CLASS CONSTANTS ==========
	static const double AFFY_ALPHA_MAX;// = 0.49;

private: // ========== INSTANCE DATA ==========
	double _maxErrorRate;
};


// ========== GenotypeErrorModel inline IMPLEMENTATION ==========
inline Matrix<double> const& GenotypeErrorModel::getErrorModelForSnp(int snp) const {
	return _matrices[snp];
}

inline double GenotypeErrorModel::getProbability(int snp, int O, int G) const {
	return _matrices[snp][O][G];
}

inline int GenotypeErrorModel::getNbSnps() const {
	return SIZE(getMatrix());
}

inline vector< Matrix<double> >& GenotypeErrorModel::getMatrix() {
	return _matrices;
}

inline vector< Matrix<double> > const& GenotypeErrorModel::getMatrix() const {
	return _matrices;
}


// ========== SimpleGenotypeErrorModel inline IMPLEMENTATION ==========
inline double SimpleGenotypeErrorModel::getMaxErrorRate() const {
	return _maxErrorRate;
}

// ========== AffymetrixGenotypeErrorModel inline IMPLEMENTATION ==========
inline double AffymetrixGenotypeErrorModel::getMaxErrorRate() const {
	return _maxErrorRate;
}

#endif /* GENOTYPE_ERROR_MODEL_HPP_ */
