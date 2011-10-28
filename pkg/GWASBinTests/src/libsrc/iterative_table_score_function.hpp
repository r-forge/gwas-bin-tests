/*
 * iterative_table_score_function.hpp
 *
 * Define an interface for Score Functions that operate on successive univariate or divariate contingency
 * tables.
 *
 *  Created on: Sep 21, 2011
 */

#ifndef ITERATIVE_TABLE_SCORE_FUNCTION_HPP_
#define ITERATIVE_TABLE_SCORE_FUNCTION_HPP_

#include "dataset.hpp"
#include "univariate_table.hpp"
#include "divariate_table.hpp"

// see definition below
class TableData;

//class IterativeTableScoresPvaluesComputable;

/**
 * This is an abstract class that defines an interface for score functions
 * and provide some convenient functions
 */
class IterativeTableScoreFunction {
public:
	virtual ~IterativeTableScoreFunction() {}
public: // ========== ABSTRACT INTERFACE ==========
	/**
	 * update the score statistic calculation for this table (usually for sums or products)
	 *
	 * @param data the tables for one iteration of all the tables of the dataset
	 */
	virtual void update(const TableData& data) = 0;
	/**
	 * finalize the score statistic computation for the whole dataset
	 *
	 * It is called when all the iterations are done
	 *
	 * @return the score statistic for the dataset
	 *
	 */
	virtual double compute() = 0;
public: // ========== INSTANCE METHODS ==========
	//double computeScore(PvaluesEstimator& estimator);
	//double computePvalue(int nb_permutations, PvaluesEstimator& estimator);

public: // ========== CONVENIENCE METHODS: not optimal at all ==========
	double computeScore(const Dataset& dataset);
	double computePvalue(int nb_permutations, const Dataset& dataset, int seed = -1);

	//virtual PvaluesEstimator* createDefaultEstimator(const Dataset& dataset, int seed = -1) const;
};


/**
 * class that defines the content of one iteration on all the contigency univariate and divariate tables
 * of a dataset.
 */
class TableData {
public: // ========== LIFECYCLE ==========
	TableData(int nb_phenotypes);
public: // ========== INSTANCE DATA ==========
	DivariateTable div_table;
	UnivariateTable univ_table;
	/**
	 * the index of the current_snp: N.B: the divariate table is computed on SNPs (current_snp, current_snp+1)
	 * and the univariate on current_snp
	 */
	int current_snp;
	/**
	 * the total number of snps
	 */
	int nb_snps;
};

#endif /* ITERATIVE_TABLE_SCORE_FUNCTION_HPP_ */
