/*
 * pvalues.hpp
 *
 *  Created on: Jun 24, 2010
 *      Author: kforner
 */

#ifndef PVALUES_HPP_
#define PVALUES_HPP_

#include "misc.hpp"
#include "shuffle.hpp"
#include "divariate_table.hpp"
#include "univariate_table.hpp"
#include "dataset.hpp"

class PvaluesEstimator;
/*
class DivariateUnivariateData {
public: // ========== LIFECYCLE ==========
	DivariateUnivariateData(int nb_phenotypes);
public: // ========== INSTANCE DATA ==========
	DivariateTable div_table;
	UnivariateTable univ_table;
	int index;
	int current_snp;
	int nb_snps;
};


class IterativeTableScoreFunction {
public:
	virtual ~IterativeTableScoreFunction() {}
public: // ========== ABSTRACT INTERFACE ==========
	virtual void update(const DivariateUnivariateData& data) = 0;
	virtual double compute() = 0;
public: // ========== INSTANCE METHODS ==========
	double computeScore(PvaluesEstimator& estimator);
	double computePvalue(int nb_permutations, PvaluesEstimator& estimator);

public: // ========== CONVENIENCE METHODS: not optimal ==========
	double computeScore(const Dataset& dataset);
	double computePvalue(int nb_permutations, const Dataset& dataset, int seed = -1);

	virtual PvaluesEstimator* createDefaultEstimator(const Dataset& dataset, int seed = -1) const;
};
*/
/*
class PvaluesEstimator
{
public:
	virtual ~PvaluesEstimator() {}

	// ========== INTERFACE ==========
	virtual int computePvalues(int nb_permutations,  vector<IterativeTableScoreFunction*>& score_functions,
			vector<double>& pvalues) = 0;
	virtual void computeAllscores(
			vector<IterativeTableScoreFunction*>& score_functions, vector<double>& scores) = 0;
	virtual void reset_shuffling() = 0;
	virtual void shuffle() = 0;
protected: // ========== INTERFACE TO IMPLEMENT ==========

//	virtual void shuffle(vector<UBYTE>& phenotypes) = 0;

	//	double computeScore(IterativeTableScoreFunction& score_fun, const vector<
	//			UBYTE>& phenotypes);

public:
	// ========== SHORTCUT METHODS INSTANCE METHODS ==========
	double computePvalue(int nb_permutations,
			IterativeTableScoreFunction& score_fun);
	double computeScore(IterativeTableScoreFunction& score_fun);

protected:
	// ========== INSTANCE DATA ==========
	vector<IterativeTableScoreFunction*> _score_functions_tmp;
	vector<double> _double_tmp;
};

class DivariateUnivariatePvaluesEstimator_Base: public PvaluesEstimator
{
public:
	// ========== LIFECYCLE ==========
	DivariateUnivariatePvaluesEstimator_Base(const Dataset& dataset, double min_interesting_pvalue = -1, double confidence = 0.999, double precision = -1);
	virtual ~DivariateUnivariatePvaluesEstimator_Base();

public:
	// ========== INSTANCE METHODS ==========
	virtual int computePvalues(int nb_permutations, vector<
			IterativeTableScoreFunction*>& score_functions,
			vector<double>& pvalues);
	virtual void computeAllscores(
			vector<IterativeTableScoreFunction*>& score_functions,
			vector<double>& scores);
	virtual void reset_shuffling();

protected:


	void computeAllscores_nSnps(
			vector<IterativeTableScoreFunction*>& score_functions,
			const vector<UBYTE>& encoded_phenotypes, vector<double>& scores);
	void computeAllscores_oneSnp(
			vector<IterativeTableScoreFunction*>& score_functions,
			const vector<UBYTE>& phenotypes, vector<double>& scores);

	vector<UBYTE>& getPhenotypesToShuffle();

protected:
	// ========== CONSTANTS ==========
	static const int LOOP_PERMUTATION_NUMBER = 100;
	static const double PRECISION_CONFIDENCE = 1 - 0.01;

protected:
	// ========== INSTANCE DATA ==========
	DivariateUnivariateData _table_data;
	const Dataset& _dataset;
	DivariateOptimizedDataSet* _optimized_dataset_ptr;
	vector<UBYTE> _phenotypes_to_shuffle;
	double _min_pvalue;
	double _confidence;
	double _precision;

private: // ========== FORBIDDEN METHODS ==========
	DivariateUnivariatePvaluesEstimator_Base(DivariateUnivariatePvaluesEstimator_Base const& source);
	DivariateUnivariatePvaluesEstimator_Base& operator = (DivariateUnivariatePvaluesEstimator_Base const& source);
};




template<class SHUFFLER>
class DivariateUnivariatePvaluesEstimator : public DivariateUnivariatePvaluesEstimator_Base {
public:
	DivariateUnivariatePvaluesEstimator(const Dataset& dataset, SHUFFLER& shuffler, double min_interesting_pvalue = -1, double confidence = 0.999, double precision = -1)
		: DivariateUnivariatePvaluesEstimator_Base(dataset, min_interesting_pvalue, confidence, precision), _shuffler(shuffler)  {}
protected:
	virtual void shuffle() {
		vector<UBYTE>& phenotypes = getPhenotypesToShuffle();
		_shuffler.shuffle( &phenotypes[0], SIZE(phenotypes));
	}
private:
	SHUFFLER _shuffler;
};


inline PvaluesEstimator* IterativeTableScoreFunction::createDefaultEstimator(const Dataset& dataset, int seed) const
{
	StdShuffler shuffler;
	if ( seed != -1)
		shuffler.reseed(seed);
	return new DivariateUnivariatePvaluesEstimator<StdShuffler>(dataset, shuffler);
}

//inline double IterativeTableScoreFunction::computeScore(
//		PvaluesEstimator& estimator, const vector<UBYTE>& phenotypes)
//{
//	return estimator.computeScore(*this, phenotypes);
//}

inline double IterativeTableScoreFunction::computePvalue(int nb_permutations, PvaluesEstimator& estimator)
{
	return estimator.computePvalue(nb_permutations, *this);
}

inline double IterativeTableScoreFunction::computeScore(PvaluesEstimator& estimator)
{
	return estimator.computeScore(*this);
}

inline void DivariateUnivariatePvaluesEstimator_Base::computeAllscores(vector<
		IterativeTableScoreFunction*>& score_functions,
		vector<double>& scores)
{
	if (_dataset.nb_snps() == 1)
		computeAllscores_oneSnp(score_functions, _phenotypes_to_shuffle, scores);
	else
		computeAllscores_nSnps(score_functions, _phenotypes_to_shuffle, scores);
}

inline vector<UBYTE>& DivariateUnivariatePvaluesEstimator_Base::getPhenotypesToShuffle()
{
	return _phenotypes_to_shuffle;
}

inline void DivariateUnivariatePvaluesEstimator_Base::reset_shuffling()
{
	_phenotypes_to_shuffle = _dataset.nb_snps() > 1
			? _optimized_dataset_ptr->encodedPhenotypes(): _dataset.phenotypes();

}*/

#endif /* PVALUES_HPP_ */
