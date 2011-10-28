/*
 * iterative_table_score_function.cpp
 *
 *  Created on: Sep 21, 2011
 */

#include "iterative_table_score_function.hpp"
#include "pvalues_computer.hpp"
#include "iterative_table_scores_pvalues_computable.hpp"

TableData::TableData(int nb_phenotypes)
	: div_table(nb_phenotypes), univ_table(nb_phenotypes), current_snp(0)
{
}

double IterativeTableScoreFunction::computePvalue(int nb_permutations, const Dataset& dataset, int seed)
{
	PvaluesComputer pvc;
	StdShuffler shuffler(seed, 0);
	IterativeTableScoresPvaluesComputable comp(*this, dataset, shuffler);
	pvc.computePvalues(nb_permutations, comp);
	return pvc.get_pvalues().front();
}

double IterativeTableScoreFunction::computeScore(const Dataset& dataset)
{
	StdShuffler shuffler(0, 0);
	IterativeTableScoresPvaluesComputable comp(*this, dataset, shuffler);
	vector<double> scores;

	comp.compute( scores );
	return scores.front();
}
