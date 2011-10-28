/*
 * group_iterative_table_scores_pvalues_computable.hpp
 *
 * An extension of <IterativeTableScoresPvaluesComputable> that handles groups
 * i.e the statistics and permutations will be done within groups and
 * the sum of groups score will be transparently computed.
 *
 * It will internally use IterativeTableScoresPvaluesComputable instances for each group.
 * The shuffler will be shared by each IterativeTableScoresPvaluesComputable instance.
 *
 *  Created on: Sep 22, 2011
 *
 */

#ifndef GROUP_ITERATIVE_TABLE_SCORES_PVALUES_COMPUTABLE_HPP_
#define GROUP_ITERATIVE_TABLE_SCORES_PVALUES_COMPUTABLE_HPP_

#include "iterative_table_scores_pvalues_computable.hpp"
#include "group_manager.hpp"

class GroupIterativeTableScoresPvaluesComputable : public PvaluesComputable {
public: // ==== IMPLEMENTATION of PvaluesComputable INTERFACE=====
	virtual void compute(vector<double>& scores);
	virtual void shuffle();
public: // ===== LIFECYCLE =====
	GroupIterativeTableScoresPvaluesComputable(vector<IterativeTableScoreFunction*> score_functions,
			const Dataset& dataset, const vector<int>& groups, StdShuffler& shuffler);
	GroupIterativeTableScoresPvaluesComputable(IterativeTableScoreFunction& score_function,
			const Dataset& dataset, const vector<int>& groups, StdShuffler& shuffler);

	virtual ~GroupIterativeTableScoresPvaluesComputable();

protected:
	void init(vector<IterativeTableScoreFunction*> score_functions, const Dataset& dataset, const vector<int>& groups, StdShuffler& shuffler);

private: // ===== INSTANCE DATA =====
	GroupManager _group_manager;
	vector<Dataset*> _subgroup_datasets;
	vector<IterativeTableScoresPvaluesComputable*> _subgroups_computable;
	vector<double> _subgroup_scores;
};

#endif /* GROUP_ITERATIVE_TABLE_SCORES_PVALUES_COMPUTABLE_HPP_ */
