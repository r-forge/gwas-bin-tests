/*
 * pvalues_computer.cpp
 *
 *  Created on: Sep 20, 2011
 */

#include "pvalues_computer.hpp"
#include "binomial.hpp"

INT64 PvaluesComputer::computePvalues(INT64 nb_permutations, PvaluesComputable & computable)
{
	const bool use_minpvalue_heuristic = _min_pvalue > 0;
	const bool use_max_relative_error_heuristic = _max_relative_error > 0;
	const bool use_heuristic = use_minpvalue_heuristic || use_max_relative_error_heuristic;
	const int freq = _heuristic_loop_frequency;

	//computable.init();
	computable.compute( _observed_scores );

	const int nb_scores = SIZE(_observed_scores);

	vector<int> nb_better(nb_scores, 0);
	vector<double> scores;

	INT64 k = 0;
	for (k = 0; k < nb_permutations; ++k)
	{
		// heuristics: check if we can stop the calculus
		if (use_heuristic && k > 0 && (k % freq == 0)) {
			if (use_minpvalue_heuristic && check_minpvalue_heuristic(_min_pvalue, k, nb_better, _confidence) )
				break;
			if ( use_max_relative_error_heuristic && check_max_relative_error_heuristic(_max_relative_error, k, nb_better, _confidence) )
				break;
		}

		computable.shuffle();
		scores.resize(0);
		computable.compute( scores );

		// update nb_better
		for (int i = 0; i < nb_scores; ++i)
			if ( scores[i] >= _observed_scores[i] )
				nb_better[i]++;
	}

	// compute the pvalues: pv = (#{scores >= observed}+1) / (nb_permutations + 1)
	_pvalues.resize(0);
	for (int i = 0; i < nb_scores; ++i)
		_pvalues.push_back( (nb_better[i] + 1) / (k + 1.0) );

	return k;
}

bool PvaluesComputer::check_minpvalue_heuristic(double min_pvalue, INT64 nb_permutations, const vector<int> & nb_better, double confidence)
{
	const int nb = SIZE( nb_better );
	bool all_above = true;
	for (int i = 0; i < nb; ++i) {
		if ( ! Binomial::is_pvalue_above(nb_better[i]/double(nb_permutations), min_pvalue, nb_permutations, confidence) ) {
			all_above = false;
			break;
		}
	}

	return all_above;
}

bool PvaluesComputer::check_max_relative_error_heuristic(double max_relative_error, INT64 nb_permutations, const vector<int> & nb_better, double confidence)
{
	const int nb = SIZE( nb_better );
	bool all_within_precision = true;
	double error = 0;
	for (int i = 0; i < nb; ++i) {
		error = Binomial::max_relative_error(nb_better[i]/double(nb_permutations), nb_permutations, confidence);
		if ( error > max_relative_error ) {
			all_within_precision = false;
			break;
		}
	}
	return all_within_precision;
}

std::ostream & PvaluesComputer::printOn(std::ostream & o) const
{
	o << "PvaluesComputer(";
	o << "min_pvalue=" << get_min_pvalue();
	o << ", max_relative_error=" << get_max_relative_error();
	o << ", confidence=" << get_confidence();
	o << ", heuristic_loop_frequency=" << get_heuristic_loop_frequency();
	o <<")\n";

	if ( ! _observed_scores.empty() ) {
		o << "Last Observed scores=";
		for (int i = 0; i < SIZE(_observed_scores); ++i)
			o << _observed_scores[i] << ", ";
		o << endl;
	} else
		o << "No Observed Scores calculated yet\n";

	if ( ! _pvalues.empty() ) {
		o << "Last pvalues=";
		for (int i = 0; i < SIZE(_pvalues); ++i)
			o << _pvalues[i] << ", ";
		o << endl;
	} else
		o << "No pvalues calculated yet\n";
	o << endl;

	return o;
}
