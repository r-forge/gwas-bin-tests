/*
 * univariate_score_function.cpp
 *
 *  Created on: Sep 21, 2011
 *      Author: m160508
 */

#include "univariate_score_function.hpp"
#include "genotype_error_model.hpp"
UnivariateScoreFunction::UnivariateScoreFunction(int nb_phenotypes, double regCoeff, GenotypeErrorModel const& error_models)
 : _univ_llh(nb_phenotypes, regCoeff),
   _error_models( error_models ),
   _current_dataset_likelihood(0)
{
}

void UnivariateScoreFunction::update(const TableData& data) {
	double score = _univ_llh.logLikelihood(data.univ_table, _error_models.getErrorModelForSnp(data.current_snp), true);
	_current_dataset_likelihood += score; // partial sum
}

double UnivariateScoreFunction::compute() {
	double dataset_likelihood = _current_dataset_likelihood;
	_current_dataset_likelihood = 0;
	return dataset_likelihood;
}
