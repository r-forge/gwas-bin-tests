/*
 * univariate_score_function.hpp
 *
 * Implements a IterativeTableScoreFunction class for the UnivariateLikelihood statistic.
 *
 *  Created on: Sep 21, 2011

 */

#ifndef UNIVARIATE_SCORE_FUNCTION_HPP_
#define UNIVARIATE_SCORE_FUNCTION_HPP_

#include "iterative_table_score_function.hpp"
#include "univariate_likelihood.hpp"

class UnivariateScoreFunction : public IterativeTableScoreFunction {
public: // ========== IMPLEMENTATION of IterativeTableScoreFunction INTERFACE ==========
	virtual void update(const TableData& data);
	virtual double compute();

public: // ========== LIFECYCLE ==========
	UnivariateScoreFunction(int nb_phenotypes, double regCoeff, GenotypeErrorModel const& error_models);
	virtual ~UnivariateScoreFunction() {}

private: // ========== INSTANCE DATA ==========
	UnivariateLikelihood _univ_llh;
	GenotypeErrorModel const& _error_models;
	double _current_dataset_likelihood;
};


#endif /* UNIVARIATE_SCORE_FUNCTION_HPP_ */
