/*
 * iterative_table_scores_pvalues_computable.cpp
 *
 *  Created on: Sep 21, 2011
 */

#include "iterative_table_scores_pvalues_computable.hpp"

IterativeTableScoresPvaluesComputable::IterativeTableScoresPvaluesComputable(IterativeTableScoreFunction& score_function,
		const Dataset& dataset, StdShuffler& shuffler)
		: PhenotypesShufflingPvaluesComputable(dataset.phenotypes(), shuffler),
		  _score_functions(1, &score_function),
		  _table_data(dataset.nb_phenotypes()),
		  _dataset(dataset),
		  _optimized_dataset_ptr(0)
{
	_table_data.nb_snps = dataset.nb_snps();

	// if more than one SNP, we use an optimized dataset that enables
	// quick divariate table computation. It uses encoded phenotypes too.
	if (dataset.nb_snps() > 1) {
		_optimized_dataset_ptr = new DivariateOptimizedDataSet(dataset);
		set_phenotypes( _optimized_dataset_ptr->encodedPhenotypes() );
	}
}

IterativeTableScoresPvaluesComputable::IterativeTableScoresPvaluesComputable(vector<IterativeTableScoreFunction*> score_functions,
		const Dataset& dataset, StdShuffler& shuffler)
		: PhenotypesShufflingPvaluesComputable(dataset.phenotypes(), shuffler),
		  _score_functions(score_functions),
		  _table_data(dataset.nb_phenotypes()),
		  _dataset(dataset),
		  _optimized_dataset_ptr(0)
{
	_table_data.nb_snps = dataset.nb_snps();

	// if more than one SNP, we use an optimized dataset that enables
	// quick divariate table computation. It uses encoded phenotypes too.
	if (dataset.nb_snps() > 1) {
		_optimized_dataset_ptr = new DivariateOptimizedDataSet(dataset);
		set_phenotypes( _optimized_dataset_ptr->encodedPhenotypes() );
	}
}

IterativeTableScoresPvaluesComputable::~IterativeTableScoresPvaluesComputable()
{
	if ( _optimized_dataset_ptr )
		delete _optimized_dataset_ptr;
}

void IterativeTableScoresPvaluesComputable::reset_shuffling()
{
	set_phenotypes( _optimized_dataset_ptr ? _optimized_dataset_ptr->encodedPhenotypes()
			: _dataset.phenotypes() );
}

void IterativeTableScoresPvaluesComputable::compute_for_phenotypes(const vector<UBYTE> & phenotypes, vector<double> & scores)
{
	if (_dataset.nb_snps() == 1)
		compute_oneSnp(phenotypes, scores);
	else
		compute_manySnps(phenotypes, scores);
}

void IterativeTableScoresPvaluesComputable::compute_oneSnp(const vector<UBYTE> & phenotypes, vector<double> & scores)
{
	assert( _table_data.nb_snps == _dataset.nb_snps() && _dataset.nb_snps() == 1);

	const int nb_scores = SIZE(_score_functions);

	scores.resize(0);

	_table_data.current_snp = 0;
	_table_data.univ_table.fill(_dataset.genotypes()[0], phenotypes);
	for (int i = 0; i < nb_scores; ++i) {
		_score_functions[i]->update(_table_data);
		scores.push_back( _score_functions[i]->compute() );
	}
}

void IterativeTableScoresPvaluesComputable::compute_manySnps(const vector<UBYTE> & phenotypes, vector<double> & scores)
{
	assert( _dataset.nb_snps() > 1 && _optimized_dataset_ptr != 0);
	assert( _dataset.nb_snps() == _optimized_dataset_ptr->nb_snps());
	assert( _table_data.nb_snps == _dataset.nb_snps());

	const int nb_scores = SIZE(_score_functions);
	const int nb_snps_pairs = _optimized_dataset_ptr->nb_snp_pairs();
	scores.resize(0);

	DivariateTable& div_table = _table_data.div_table;
	UnivariateTable& univ_table = _table_data.univ_table;

	for (int snp = 0; snp < nb_snps_pairs; ++snp)
	{
		_table_data.current_snp = snp;
		div_table.fill(_optimized_dataset_ptr->encodedGenotypePairs()[snp], phenotypes);
		univ_table.fillForFirstSnp(div_table);
		for (int i = 0; i < nb_scores; ++i)
			_score_functions[i]->update(_table_data);
	}

	// last snp for scores based on univariate tables
	_table_data.current_snp = nb_snps_pairs;
	univ_table.fillForSecondSnp(div_table);
	for (int i = 0; i < nb_scores; ++i)
		_score_functions[i]->update(_table_data);

	for (int i = 0; i < nb_scores; ++i)
		scores.push_back( _score_functions[i]->compute() );
}
