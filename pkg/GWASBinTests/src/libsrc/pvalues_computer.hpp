/*
 * pvalues_computer.hpp
 *
 * This class implements a generic permutation-based pvalues calculation
 * algorithm supporting two optional heuristics to speed-up. It does not implement any score
 * nor any shuffling by itself, just the general pvalues calculation and the associated heuristics.
 *
 * Heuristics:
 *
 * 1) The min_pvalue "minimum interesting pvalue" heuristic
 *
 * Principle: pvalues above the min_pvalue threshold (with confidence confidence)
 * based on their intermediate value may be computed with less permuations
 * 		The goal of this strategy is to avoid speeding permutations on getting a good accuracy on p-values that
 * 		are of no interest. Let's say for instance we are computing p-values for 1000,000 SNPs, we for sure
 * 		are not interested in p-values > 0.1 .
 * 		The estimated p-value follows a binomial distribution, so that we can compute its lower bound for a
* 		given probability, controlled by the parameter "confidence".
* 		So if the lower bound of the p-value is greater than "min_pvalue" we stop the computation,
* 		and report the current estimation of the p-value along with the actual number of permutations used to compute it.
*
 * 2) the "max_relative_error" heuristic
 *
* 		Basically the principle of this heuristic is that the lower the p-value, the higher the number of permutations
* 		it needs to get a good accuracy, and most often we are greatly interested in those small p-values.
* 		What precision do we need ? We could set a maximum value on the radius of the confidence interval
* 		of the estimated p-value, e.g. 10^6. But what if the real p-value is 0.01 ?
* 		Then this amount of precision is a waste. And what if the real p-value is \10^-7 ?
* 		Then knowing that its confidence interval radius is < 10^-6	 is of not great use.
*
* 		So instead the idea is to control the relative error on the approximated p-value
		Given the parameters max_relative_error confidence
* 		the computation will stop when the relative error at a given confidence level will be lower than "max_relative_error"
 *
 *  Created on: Sep 20, 2011
 *
 */

#ifndef PVALUES_COMPUTER_HPP_
#define PVALUES_COMPUTER_HPP_

#include "misc.hpp"
#include "pvalues_computable.hpp"

class PvaluesComputer {
public: // ========== LIFECYCLE ==========
	/**
	 * build an instance with the heuristics disabled and a default confidence of 99.9%
	 *
	 * See set_min_pvalue and set_max_relative_error to activate the heuristics
	 *
	 */
	PvaluesComputer() :
			_min_pvalue(-1), _confidence(0.999), _max_relative_error(-1), _heuristic_loop_frequency(100) {}
public:
	// ========== ACCESSORS ==========
	/* will return the last computed pvalues. If <computePvalues> has not yet
	 * been executed it will return an empty vector
	 */
	const vector<double>& get_pvalues() const { return _pvalues; }
	/* will return the last computed scores. If <computePvalues> has not yet
	 * been executed it will return an empty vector
	 */
	const vector<double>& get_observed_scores() const { return _observed_scores; }
	double get_min_pvalue() const { return _min_pvalue; }
	double get_confidence() const { return _confidence; }
	double get_max_relative_error() const { return _max_relative_error; }
	int get_heuristic_loop_frequency() const { return _heuristic_loop_frequency; }

public:
	// ========== SETTERS ==========
	/**
	 * Set the min_pvalue for the "minimum interesting pvalue" heuristic.
	 */
	void set_min_pvalue(double min_pvalue)  {  _min_pvalue = min_pvalue; }

	/**
	 * Set the max_relative_error for the heuristic (see above)
	 * a value of <= 0 disables the heuristic.
	 */
	void set_max_relative_error(double error)  {  _max_relative_error = error ; }

	/**
	 * Set the confidence for the confidence intervals used by both heuristics
	 * The default value is 0.999 (99.9%)
	 */
	void set_confidence(double confidence)  {  _confidence = confidence; }

	/** set the heuristic loop frequency, i.e the number of permutations afer which
	 * the heuristics will be checked in order to see if the calculations can stop
	 *
	 * The default is 100
	 *
	 */
	void set_heuristic_loop_frequency(int freq)  { _heuristic_loop_frequency = freq; }



public:
	// ========== INSTANCE METHODS ==========

	/** compute the pvalues for the computable
	 * @param nb_permutations the maximum number of permutations to use
	 * @param computable to compute the scores and manage the shuffling
	 *
	 * @return the maximum number of permutations done for all the pvalues
	 */
	INT64 computePvalues(INT64 nb_permutations, PvaluesComputable& computable);

	std::ostream& printOn(std::ostream& o) const;
	friend std::ostream& operator<<(std::ostream& o, PvaluesComputer const& c) {
		return c.printOn(o);
	}

public:
	// ========== CLASS METHODS ==========
	/**
	 * implements the "minimum interesting pvalue" heuristic
	 *
	 * @param min_pvalue the minimum interesting pvalue threshold
	 * @param nb_permutations the number of permutations done so far
	 * @param nb_better a vector of the number of times a better score than the observed has been found
	 * @param confidence the confidence to use for the lower bound estimation
	 *
	 * @return true iff all the pvalues are above the threshold (and therefore are not interesting
	 * with confidence "confidence"
	 */
	bool check_minpvalue_heuristic(double min_pvalue, INT64 nb_permutations, const vector<int>& nb_better, double confidence);

	/**
	 * implements the "max_relative_error" heuristic
	 *
	 * @param max_relative_error the minimum precision threshold
	 * @param nb_permutations the number of permutations done so far
	 * @param nb_better a vector of the number of times a better score than the observed has been found
	 * @param confidence the confidence to use for the lower bound estimation
	 *
	 * @return true iff all the pvalues have a precision above the threshold (and therefore are precise enough
	 * with confidence "confidence"
	 */
	bool check_max_relative_error_heuristic(double max_relative_error,  INT64 nb_permutations, const vector<int>& nb_better, double confidence);

protected:
	// ========== INSTANCE DATA ==========
	vector<double> _pvalues;
	vector<double> _observed_scores;
	double _min_pvalue;
	double _confidence;
	double _max_relative_error;
	int _heuristic_loop_frequency;
};


#endif /* PVALUES_COMPUTER_HPP_ */
