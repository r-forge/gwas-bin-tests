/*
 * divariate_score_function.cpp
 *
 *  Created on: Sep 21, 2011
 *      Author: m160508
 */

#include "divariate_score_function.hpp"
#include "genotype_error_model.hpp"

DivariateScoreFunction::DivariateScoreFunction(int nb_phenotypes, double regCoeff,
		GenotypeErrorModel const& error_models)
	: _div_llh(nb_phenotypes, regCoeff),
	  _error_models( error_models ),
	  _current_dataset_likelihood(0)
{}

// ========== IMPLEMENTATION of IterativeTableScoreFunction INTERFACE ==========

void DivariateScoreFunction::update(const TableData& data)
{
	if ( data.current_snp + 1 >= data.nb_snps)
		return;
	double l = _div_llh.logLikelihood(data.div_table, _error_models.getErrorModelForSnp(
			data.current_snp), _error_models.getErrorModelForSnp(
					data.current_snp + 1), true);

	_current_dataset_likelihood += l;
}

double DivariateScoreFunction::compute()
{
	double dataset_likelihood = _current_dataset_likelihood;
	_current_dataset_likelihood = 0;
	return dataset_likelihood;
}

