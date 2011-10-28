/*
 * genotype_error_model.cpp
 *
 *  Created on: Jun 18, 2010
 *      Author: kforner
 */

#include "genotype_error_model.hpp"
#include <cmath>




GenotypeErrorModel::GenotypeErrorModel(const vector<vector<UBYTE> > & data)
	: _matrices(SIZE(data), Matrix<double>(NB_GENOTYPES_WITH_ZZ, NB_GENOTYPES_NO_ZZ))
{
}

const double SimpleGenotypeErrorModel::DEFAULT_ERROR_RATE = 0.05;

SimpleGenotypeErrorModel::SimpleGenotypeErrorModel(const vector<vector<UBYTE> > & data,
		double max_error_rate) :
	GenotypeErrorModel(data), _maxErrorRate(max_error_rate) {
	const int nb_snps = SIZE(data);
	const int nb_samples = SIZE(data[0]);
	double alpha = max_error_rate;
	int counts[NB_GENOTYPES_WITH_ZZ];

	for (int i = 0; i < nb_snps; ++i)
	{
		count_genotypes(data[i], counts);
		double beta = counts[MISSING] / (double) nb_samples;
		fillGenotypeLikelihoodMatrixForSimpleErrorModel(alpha, beta, getMatrix()[i]);
	}
}



void SimpleGenotypeErrorModel::fillGenotypeLikelihoodMatrixForSimpleErrorModel(
		double alpha, double beta, Matrix<double>& P)
{
	double same = (1 - beta) * (1 - alpha);
	double error = (1 - beta) * alpha / 2.0;

	/*  line by line */

	/*  first line */
	P[HOMOZYGOTE1][HOMOZYGOTE1] = same;
	P[HOMOZYGOTE1][HETEROZYGOTE] = error;
	P[HOMOZYGOTE1][HOMOZYGOTE2] = error;

	/*  second line */
	P[HETEROZYGOTE][HOMOZYGOTE1] = error;
	P[HETEROZYGOTE][HETEROZYGOTE] = same;
	P[HETEROZYGOTE][HOMOZYGOTE2] = error;

	/*  3rd line */
	P[HOMOZYGOTE2][HOMOZYGOTE1] = error;
	P[HOMOZYGOTE2][HETEROZYGOTE] = error;
	P[HOMOZYGOTE2][HOMOZYGOTE2] = same;

	/*  last line */
	P[MISSING][HOMOZYGOTE1] = beta;
	P[MISSING][HETEROZYGOTE] = beta;
	P[MISSING][HOMOZYGOTE2] = beta;
}


std::ostream & operator <<(std::ostream & o, const GenotypeErrorModel & m)
{
	o << "GenotypeErrorModel for " << m.getNbSnps() << endl;
	o << "Model for first SNP:\n" << m.getErrorModelForSnp(0);
	return o;
}


void GenotypeErrorModel::count_genotypes(const vector<UBYTE> & genotypes, int counts[NB_GENOTYPES_WITH_ZZ])
{
	for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
		counts[i] = 0;
	foreach(UBYTE g, genotypes)
		counts[g]++;

}

const double AffymetrixGenotypeErrorModel::AFFY_ALPHA_MAX = 0.49;

AffymetrixGenotypeErrorModel::AffymetrixGenotypeErrorModel(const vector<vector<UBYTE> > & data,
		double max_error_rate) :
	GenotypeErrorModel(data), _maxErrorRate(max_error_rate) {

	const int nb_snps = SIZE(data);
	double alpha = max_error_rate;
	int counts[NB_GENOTYPES_WITH_ZZ];

	for (int i = 0; i < nb_snps; ++i)
	{
		// TODO: check the value of alpha, different from nnbc-Code
		alpha = max_error_rate;
		count_genotypes(data[i], counts);
		double beta = solveAffymetrixMissingRate(alpha, counts);
		if (beta < 0)
			fixAffymetrixRatesWhenNoSolution(counts, alpha, beta);
		else
			alpha = max_error_rate;
		fillGenotypeLikelihoodMatrixForAffymetrixModel(alpha, beta, getMatrix()[i]);
	}
}

void AffymetrixGenotypeErrorModel::fillGenotypeLikelihoodMatrixForAffymetrixModel(double alpha, double beta,
		Matrix<double>& P)
{

	double em = (1 - alpha) * (1 - beta);
	double cem = alpha * (1 - beta);
	double cem2 = (1 - beta - beta) * alpha;
	double e2m2 = (1 - alpha - alpha) * (1 - beta - beta);
	double m2 = beta + beta;

	/* line by line */

	/*// first line*/
	P[HOMOZYGOTE1][HOMOZYGOTE1] = em;
	P[HOMOZYGOTE1][HETEROZYGOTE] = cem2;
	P[HOMOZYGOTE1][HOMOZYGOTE2] = 0;

	/*// second line*/
	P[HETEROZYGOTE][HOMOZYGOTE1] = cem;
	P[HETEROZYGOTE][HETEROZYGOTE] = e2m2;
	P[HETEROZYGOTE][HOMOZYGOTE2] = cem;

	/*// 3rd line*/
	P[HOMOZYGOTE2][HOMOZYGOTE1] = 0;
	P[HOMOZYGOTE2][HETEROZYGOTE] = cem2;
	P[HOMOZYGOTE2][HOMOZYGOTE2] = em;

	/*// last line*/
	P[MISSING][HOMOZYGOTE1] = beta;
	P[MISSING][HETEROZYGOTE] = m2;
	P[MISSING][HOMOZYGOTE2] = beta;
}


void AffymetrixGenotypeErrorModel::fixAffymetrixRatesWhenNoSolution(int counts[NB_GENOTYPES_WITH_ZZ],
		double& alpha_out, double& beta_out, bool debug)
{

	double size = counts[HOMOZYGOTE1]+counts[HETEROZYGOTE]+counts[HOMOZYGOTE2]+counts[MISSING];
	double pmo = counts[MISSING] / size;
	double pho = counts[HETEROZYGOTE] / size;
	beta_out = pmo;
	/*And the adjusted alpha is :*/
	if (pmo < 1.0)
		alpha_out = pho / (1 - pmo);
	else
		/*If there are only missing values...
		 Usually, such a feature is eliminated before*/
		alpha_out = 1.0;

	if (alpha_out > AFFY_ALPHA_MAX) {
		alpha_out = AFFY_ALPHA_MAX; /* avoid zeroes in the table */
		if (debug)
		{
			fprintf(stderr,"fixAffymetrixRatesWhenNoSolution():setting alpha to %f because it was too high\n", alpha_out);

		}
	}

	if (debug > 1)
	{
		fprintf(
				stderr,
				"fixAffymetrixRatesWhenNoSolution():s probably due to a SNP that should have been filtered\n");
		fprintf(stderr, "hom1=%i hom2=%i het=%i ZZ=%i total=%i\n",
				(int) counts[HOMOZYGOTE1], (int) counts[HOMOZYGOTE2],
				(int) counts[HETEROZYGOTE], (int) counts[MISSING],
				(int) size);
		fprintf(stderr, "ZZ rate is %f\n", (float) (counts[MISSING]
				/ (size + 1.0)));
		fprintf(stderr, "Computed alpha=%f beta=%f\n", (float) alpha_out,
				(float) beta_out);
	}

}

double AffymetrixGenotypeErrorModel::solveAffymetrixMissingRate(double alpha, int counts[NB_GENOTYPES_WITH_ZZ])
{
	int i;
	double size, pho, pmo, aa, bb, cc, dd, roots[2], r, z;
	double EPS = 1e-6;
	if (counts[MISSING] <= 0)
		return 0;

	/* the unknown value to solve (i.e the beta) */

	size = counts[HOMOZYGOTE1]+counts[HETEROZYGOTE]+counts[HOMOZYGOTE2]+counts[MISSING];
	/*Proportion of heterozygotous Noted P(Ojb = Aa) in the article*/
	pho = counts[HETEROZYGOTE] / size;

	/*Proportion of missing Noted P(Ojb = 0) in the article*/
	pmo = counts[MISSING] / size;

	/*Replacing z by pmo/beta - 1 in Eq 17 leads to the following second degree equation
	 aa x beta^2 + bb + cc =0*/
	aa = 2 * (1 - 3 * alpha);
	bb = -pho + pmo * (-2 + 5 * alpha) - 1 + 4 * alpha;
	cc = pmo * (1 - 3 * alpha);

	/*Let's us solve it.Its discriminant is*/
	dd = bb * bb - 4 * aa * cc;

	if (dd < 0.0) /* no root */
		return -1.0;

	dd = sqrt(dd);
	/*The two roots are*/
	roots[0] = (-bb - dd) / (aa + aa);
	roots[1] = (-bb + dd) / (aa + aa);

	for (i = 0; i < 2; ++i)
	{
		r = roots[i];
		if (r > 0 && r <= (0.5 + EPS))
		{
			z = pmo / r - 1;
			if (z >= -EPS && z <= 1 + EPS)
			{
				return r;
			}
		}
	}

	return -2.0;
}


