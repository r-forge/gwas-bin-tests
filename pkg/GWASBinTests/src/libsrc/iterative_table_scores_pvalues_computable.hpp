/*
 * iterative_table_scores_pvalues_computable.hpp
 *
 * Define a PvaluesComputable class that uses a list of IterativeTableScoreFunctions
 * to compute the scores. It handles the shuffling of phenotypes via the PhenotypesShufflingPvaluesComputable
 * helper class and iterate over the data set to supply the TableData to the score functions.
 *
 *  Created on: Sep 21, 2011
 */

#ifndef ITERATIVE_TABLE_SCORES_PVALUES_COMPUTABLE_HPP_
#define ITERATIVE_TABLE_SCORES_PVALUES_COMPUTABLE_HPP_

#include "iterative_table_score_function.hpp"
#include "pvalues_computable.hpp"

class IterativeTableScoresPvaluesComputable : public PhenotypesShufflingPvaluesComputable {
public: // ==== IMPLEMENTATION OF PhenotypesShufflingPvaluesComputable INTERFACE =====
	/**
	 * This is the main function, it will iterate and build all tables and call each score function on them
	 */
	virtual void compute_for_phenotypes(const vector<UBYTE>& phenotypes, vector<double>& scores);

public: // ===== LIFECYCLE =====
	IterativeTableScoresPvaluesComputable(IterativeTableScoreFunction& score_function,
				const Dataset& dataset, StdShuffler& shuffler);
	IterativeTableScoresPvaluesComputable(vector<IterativeTableScoreFunction*> score_functions,
			const Dataset& dataset, StdShuffler& shuffler);
	virtual ~IterativeTableScoresPvaluesComputable();

public: // ===== ACCESSORS =====
	/**
	 * the number of scores computable == the number of score functions
	 */
	int nb_scores() const { return SIZE(_score_functions); }

protected:
	/**
	 * implementation of compute() when there is only one SNP in the dataset
	 */
	void compute_oneSnp(const vector<UBYTE>& phenotypes, vector<double>& scores);
	/**
	 * implementation of compute() when in the general case
	 */
	void compute_manySnps(const vector<UBYTE>& phenotypes, vector<double>& scores);

	void reset_shuffling();

	// ========== INSTANCE DATA ==========
	vector<IterativeTableScoreFunction*> _score_functions;
	TableData _table_data;
	const Dataset& _dataset;
	DivariateOptimizedDataSet* _optimized_dataset_ptr;
};

#endif /* ITERATIVE_TABLE_SCORES_PVALUES_COMPUTABLE_HPP_ */
