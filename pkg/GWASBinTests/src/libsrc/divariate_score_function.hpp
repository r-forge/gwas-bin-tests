/*
 * divariate_score_function.hpp
 *
 *  Created on: Sep 21, 2011
 *      Author: m160508
 */

#ifndef DIVARIATE_SCORE_FUNCTION_HPP_
#define DIVARIATE_SCORE_FUNCTION_HPP_

#include "iterative_table_score_function.hpp"
#include "divariate_likelihood.hpp"

class DivariateScoreFunction : public IterativeTableScoreFunction {
public: // ========== IMPLEMENTATION of IterativeTableScoreFunction INTERFACE ==========
	virtual void update(const TableData& data);
	virtual double compute();

public: // ========== LIFECYCLE ==========
	DivariateScoreFunction(int nb_phenotypes, double regCoeff, GenotypeErrorModel const& error_models);
	virtual ~DivariateScoreFunction() {}

private: // ========== INSTANCE DATA ==========
	DivariateLikelihood _div_llh;
	GenotypeErrorModel const& _error_models;
	double _current_dataset_likelihood;
};

#endif /* DIVARIATE_SCORE_FUNCTION_HPP_ */
