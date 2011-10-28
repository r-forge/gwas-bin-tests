/*
 * pvalues.cpp
 *
 *  Created on: Jun 24, 2010
 *      Author: kforner
 */

#include "pvalues.hpp"
#include "dataset.hpp"

#include <boost/math/distributions/binomial.hpp>
using namespace boost::math;
/*
double IterativeTableScoreFunction::computeScore(const Dataset& dataset)
{
	PvaluesEstimator* estimator = createDefaultEstimator(dataset);
	double score = computeScore(*estimator);
	delete estimator;
	return score;
}


double IterativeTableScoreFunction::computePvalue(int nb_permutations, const Dataset& dataset, int seed)
{
	PvaluesEstimator* estimator = createDefaultEstimator(dataset, seed);
	double pvalue = computePvalue(nb_permutations, *estimator);
	delete estimator;
	return pvalue;
}

double PvaluesEstimator::computePvalue(int nb_permutations,
		IterativeTableScoreFunction& score_fun)
{
	_score_functions_tmp.resize(0);
	_score_functions_tmp.push_back(&score_fun);
	computePvalues(nb_permutations, _score_functions_tmp, _double_tmp);
	return _double_tmp[0];
}

double PvaluesEstimator::computeScore(
		IterativeTableScoreFunction& score_fun)
{
	_score_functions_tmp.resize(0);
	_score_functions_tmp.push_back(&score_fun);
	computeAllscores(_score_functions_tmp, _double_tmp);
	return _double_tmp[0];
}


DivariateUnivariatePvaluesEstimator_Base::DivariateUnivariatePvaluesEstimator_Base(
		const Dataset & dataset, double min_interesting_pvalue, double confidence,  double precision)
 : _table_data(dataset.nb_phenotypes()), _dataset(dataset), _optimized_dataset_ptr(0),
   _min_pvalue(min_interesting_pvalue), _confidence(confidence), _precision(precision)
{
	_table_data.nb_snps = dataset.nb_snps();
	if (dataset.nb_snps() > 1)
		_optimized_dataset_ptr = new DivariateOptimizedDataSet(dataset);
	reset_shuffling();
}

DivariateUnivariatePvaluesEstimator_Base::~DivariateUnivariatePvaluesEstimator_Base()
{
	if ( _optimized_dataset_ptr )
		delete _optimized_dataset_ptr;
}

//DivariateUnivariatePvaluesEstimator_Base::DivariateUnivariatePvaluesEstimator_Base(
//				DivariateUnivariatePvaluesEstimator_Base const& source)
//	: _table_data(source._table_data), _dataset(source._dataset)
//{
//	if ( source._optimized_dataset_ptr != 0)
//		_optimized_dataset_ptr = new DivariateOptimizedDataSet(*source._optimized_dataset_ptr);
//}
//
//DivariateUnivariatePvaluesEstimator_Base
//		& DivariateUnivariatePvaluesEstimator_Base::operator =(
//				DivariateUnivariatePvaluesEstimator_Base const& source)
//{
//	if ( this != &source) {
//		_table_data = source._table_data;
//		_dataset = source._dataset;
//		if ( _optimized_dataset_ptr != 0)
//			delete _optimized_dataset_ptr;
//		if ( source._optimized_dataset_ptr != 0)
//			_optimized_dataset_ptr = new DivariateOptimizedDataSet(*source._optimized_dataset_ptr);
//	}
//	return *this;
//}

int DivariateUnivariatePvaluesEstimator_Base::computePvalues(
		int nb_permutations,
		vector<IterativeTableScoreFunction*>& score_functions
		, vector<double>& pvalues)
{
	const int nb_scores = SIZE(score_functions);

	const bool use_minpvalue_heuristic = _min_pvalue > 0;
	const bool use_precision_heuristic = _precision > 0;
	const bool use_heuristic = use_minpvalue_heuristic || use_precision_heuristic;

	const double alpha = 1.0 - _confidence;

	reset_shuffling();
//	cerr << "phenotypes_to_shuffle:\n" << phenotypes_to_shuffle;

	vector<double> observed_scores;
	_table_data.index = 0;
	computeAllscores(score_functions, observed_scores);
//	cerr << "observed=" << observed_scores[0] << endl;

	assert( SIZE(observed_scores) == nb_scores);

	vector<int> nb_better(nb_scores, 0);
	vector<double> scores;
	int k;
	for (k = 0; k < nb_permutations; ++k)
	{

		// heuristic: check if we can stop the calculus
		if (use_heuristic && k > 0 && (k % LOOP_PERMUTATION_NUMBER == 0)) {

			if (use_minpvalue_heuristic) {
				// min_pvalue heuristic
				bool all_above = true;
				for (int i = 0; i < nb_scores; ++i) {
//					cerr << k << " " << nb_better[i] << " " << alpha << endl;
					double lower =
							binomial_distribution<>::find_lower_bound_on_p(k,
									nb_better[i], alpha);
//					cerr << "lower=" << lower << endl;
					if (lower < _min_pvalue) {
						all_above = false;
						break;
					}
				}
				if (all_above)
					break;

			}

			if ( use_precision_heuristic ) {
				bool all_within_precision = true;
				for (int i = 0; i < nb_scores; ++i) {
					const int successes = nb_better[i];
					double lower =
							binomial_distribution<>::find_lower_bound_on_p(k,
									successes, alpha/2.0);
					double upper =
							binomial_distribution<>::find_upper_bound_on_p(k,
									successes, alpha/2.0);
					double diameter = upper - lower;
					double estimate = successes / double(k);
//					if ( k > 5000 )
//						cerr << "k=" << k << ", estimate=" << estimate << ", diameter=" << diameter << ", ration=" << estimate / diameter << endl;
					if (estimate / diameter  < _precision) {
						all_within_precision = false;
						break;
					}
				}
				if( all_within_precision )
					break;
			}
		}

		_table_data.index++;
		shuffle();
//		std::for_each(ALL(phenotypes_to_shuffle), cerr << (_1+0) << " ");
//		cerr << endl;
		scores.resize(0);
		computeAllscores(score_functions, scores);
//		cerr << scores[0] << endl;
		for (int i = 0; i < nb_scores; ++i)
			if ( scores[i] >= observed_scores[i] )
				nb_better[i]++;


	}

	pvalues.resize(0);
	for (int i = 0; i < nb_scores; ++i)
		pvalues.push_back( (nb_better[i] + 1) / (k + 1.0) );

	return k;
}


void DivariateUnivariatePvaluesEstimator_Base::computeAllscores_nSnps(vector<
		IterativeTableScoreFunction*>& score_functions,
		const vector<UBYTE>& encoded_phenotypes, vector<double>& scores)
{
	assert( _dataset.nb_snps() > 1 && _optimized_dataset_ptr != 0);
	assert( _dataset.nb_snps() == _optimized_dataset_ptr->nb_snps());
	assert( _table_data.nb_snps == _dataset.nb_snps());
	scores.resize(0);

	DivariateTable& div_table = _table_data.div_table;
	UnivariateTable& univ_table = _table_data.univ_table;

	const int nb_scores = SIZE(score_functions);

	const int nb_snps_pairs = _optimized_dataset_ptr->nb_snp_pairs();
	for (int snp = 0; snp < nb_snps_pairs; ++snp)
	{
		_table_data.current_snp = snp;
		div_table.fill(_optimized_dataset_ptr->encodedGenotypePairs()[snp], encoded_phenotypes);

		univ_table.fillForFirstSnp(div_table);
		for (int i = 0; i < nb_scores; ++i)
			score_functions[i]->update(_table_data);
	}

	// last snp for scores based on univariate tables
	_table_data.current_snp = nb_snps_pairs;
	univ_table.fillForSecondSnp(div_table);
	for (int i = 0; i < nb_scores; ++i)
		score_functions[i]->update(_table_data);

	for (int i = 0; i < nb_scores; ++i)
		scores.push_back( score_functions[i]->compute() );
}

void DivariateUnivariatePvaluesEstimator_Base::computeAllscores_oneSnp(vector<
		IterativeTableScoreFunction*>& score_functions,
		const vector<UBYTE>& phenotypes, vector<double>& scores)
{
	assert( _dataset.nb_snps() == 1);
	assert( _table_data.nb_snps == _dataset.nb_snps());
	scores.resize(0);
	UnivariateTable& univ_table = _table_data.univ_table;

	const int nb_scores = SIZE(score_functions);

	const int nb_snps = _dataset.nb_snps();
	for (int snp = 0; snp < nb_snps; ++snp)
	{
		_table_data.current_snp = snp;
		univ_table.fill(_dataset.genotypes()[snp], phenotypes);
		for (int i = 0; i < nb_scores; ++i)
			score_functions[i]->update(_table_data);
	}

	for (int i = 0; i < nb_scores; ++i)
		scores.push_back( score_functions[i]->compute() );
}

//
//vector<double> DivariateUnivariatePvaluesEstimator_Base::computePvalues(
//		int nb_permutations, const vector<TableScoreFunction> & score_functions) {
//
//	vector<UBYTE> phenotypes_to_shuffle = _dataset.phenotypes();
//	const int nb_scores = SIZE(score_functions);
//
//	vector<double> observed_scores;
//	for (int i = 0; i < nb_scores; ++i) {
//		observed_scores.push_back(score_functions[i].compute(_dataset))
//	}
//
//	for (int k = 0; k < nb_permutations; ++k) {
//		shuffle(phenotypes_to_shuffle);
//
//
//
//	}
//}



DivariateUnivariateData::DivariateUnivariateData(int nb_phenotypes)
	: div_table(nb_phenotypes), univ_table(nb_phenotypes), index(0)
{
}

*/
