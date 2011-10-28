/*
 * group_iterative_table_scores_pvalues_computable.cpp
 *
 *  Created on: Sep 22, 2011
 *
 */

#include "group_iterative_table_scores_pvalues_computable.hpp"



GroupIterativeTableScoresPvaluesComputable::GroupIterativeTableScoresPvaluesComputable(vector<IterativeTableScoreFunction*> score_functions, const Dataset & dataset, const vector<int> & groups, StdShuffler & shuffler)
{
	init(score_functions, dataset, groups, shuffler);
}

GroupIterativeTableScoresPvaluesComputable::GroupIterativeTableScoresPvaluesComputable(IterativeTableScoreFunction & score_function, const Dataset & dataset, const vector<int> & groups, StdShuffler & shuffler)
{
	vector<IterativeTableScoreFunction*> score_functions(1, &score_function);
	init(score_functions, dataset, groups, shuffler);
}


void GroupIterativeTableScoresPvaluesComputable::init(vector<IterativeTableScoreFunction*> score_functions, const Dataset & dataset, const vector<int> & groups, StdShuffler & shuffler)
{
	assert( _subgroup_datasets.empty() && _subgroups_computable.empty() );

	_group_manager.reset(groups);
	int nb_groups = _group_manager.nb_groups();

	for (int g = 0; g < nb_groups; ++g) {
		Dataset* subset = new Dataset(dataset, _group_manager.group_indices(g));
		_subgroup_datasets.push_back( subset );
		_subgroups_computable.push_back( new IterativeTableScoresPvaluesComputable(score_functions, *subset, shuffler) );
	}
}

GroupIterativeTableScoresPvaluesComputable::~GroupIterativeTableScoresPvaluesComputable()
{
	foreach(IterativeTableScoresPvaluesComputable* comp, _subgroups_computable)
		delete comp;
	_subgroups_computable.resize(0); // probably useless

	foreach(Dataset* set, _subgroup_datasets)
		delete set;
	_subgroup_datasets.resize(0); // probably useless
}


void GroupIterativeTableScoresPvaluesComputable::compute(vector<double> & scores)
{

	const int nb_scores = _subgroups_computable.front()->nb_scores();
	scores.resize(0); // erase vector
	scores.resize(nb_scores, 0); // init to zero

	foreach(IterativeTableScoresPvaluesComputable* comp, _subgroups_computable) {
		_subgroup_scores.resize(0);
		comp->compute(_subgroup_scores);
		assert( SIZE(_subgroup_scores) == nb_scores );
		for (int s = 0; s < nb_scores; ++s)
			scores[s] += _subgroup_scores[s];
	}
//	cerr << "scores= " << scores;
}



void GroupIterativeTableScoresPvaluesComputable::shuffle()
{
	foreach(IterativeTableScoresPvaluesComputable* comp, _subgroups_computable)
		comp->shuffle();
}

