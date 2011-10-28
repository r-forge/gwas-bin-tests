/*
 * sum_score_functions.hpp
 *
 *  Created on: Sep 21, 2011
 *      Author: m160508
 */

#ifndef SUM_SCORE_FUNCTIONS_HPP_
#define SUM_SCORE_FUNCTIONS_HPP_

#include "iterative_table_score_function.hpp"
#include "sum_scores.hpp"

class AllelicSumScoreFunction :  public IterativeTableScoreFunction {
public: // ========== IMPLEMENTATION of IterativeTableScoreFunction INTERFACE ==========
	virtual void update(const TableData& data);
	virtual double compute();
public: // ========== LIFECYCLE ==========
	AllelicSumScoreFunction(int nb_phenotypes, int nb_snps)
		:  _expected_by_snp(nb_snps, nb_phenotypes, AllelicSumScore::NB_ALLELES_NO_Z), _sum_score(0) {}
private: // ========== INSTANCE DATA ==========
	Matrix3D<double> _expected_by_snp;
	double _sum_score;
};

class GenotypicSumScoreFunction :  public IterativeTableScoreFunction {
public: // ========== IMPLEMENTATION of IterativeTableScoreFunction INTERFACE ==========
	virtual void update(const TableData& data);
	virtual double compute();
public: // ========== LIFECYCLE ==========
	GenotypicSumScoreFunction(int nb_phenotypes, int nb_snps)
		: _expected_by_snp(nb_snps, nb_phenotypes, NB_GENOTYPES_NO_ZZ), _sum_score(0) {}

private: // ========== INSTANCE DATA ==========
	Matrix3D<double> _expected_by_snp;
	double _sum_score;
};

inline void AllelicSumScoreFunction::update(const TableData & data)
{

	// TODO: the precomputation is actually useless, unless coupling it with permutations and storing
	// the expected freqs per SNP
	AllelicSumScore::precompute_allelic_expected_frequencies(data.univ_table, _expected_by_snp[data.current_snp]);
	_sum_score += AllelicSumScore::allelic_pearson_score(data.univ_table, _expected_by_snp[data.current_snp]);
}

inline double AllelicSumScoreFunction::compute()
{
	double score = _sum_score;
	_sum_score = 0;
	return score;
}

inline double GenotypicSumScoreFunction::compute()
{
	double score = _sum_score;
	_sum_score = 0;
	return score;
}


inline void GenotypicSumScoreFunction::update(const TableData & data)
{
	// TODO: the precomputation is actually useless, unless coupling it with permutations and storing
	// the expected freqs per SNP
	GenotypicSumScore::precompute_genotypic_expected_frequencies(data.univ_table, _expected_by_snp[data.current_snp]);

	_sum_score += GenotypicSumScore::genotypic_pearson_score(data.univ_table, _expected_by_snp[data.current_snp]);
}


#endif /* SUM_SCORE_FUNCTIONS_HPP_ */
