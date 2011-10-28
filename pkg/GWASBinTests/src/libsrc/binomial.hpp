/*
 * binomial.hpp
 *
 * Utility functions related to the binomial distribution (confidence intervals...)
 *
 *  Created on: Sep 20, 2011
 */

#ifndef BINOMIAL_HPP_
#define BINOMIAL_HPP_

#include "types.hpp"
#include <utility>
#include <boost/math/distributions/binomial.hpp>
using namespace boost::math;
using namespace std;

class Interval : public  pair<double, double> {
public:
	bool contains(double p) const { return p >= first && p <= second; }

	std::ostream& printOn(std::ostream& o) const { return o << "[" << first << "," << second << "]";}
	friend std::ostream& operator<<(std::ostream& o, Interval const& c) { return c.printOn(o); }
};


class Binomial {
public: // ===== STATIC/CLASS METHODS =====

	/**
	 * compute a confidence interval for either the true pvalue of for an approximated pvalue.
	 *
	 * @param the pvalue for which to compute the confidence interval
	 * @param nb_permut the number of permutation used for the approximated pvalue
	 * @param confidence the confidence of the interval
	 *
	 * @return the confidence interval as a pair of doubles
	 */
	static Interval confidence_interval(double pvalue, INT64 nb_permut,  double confidence);

	/**
	 * compute the maximum error for the approximation of the pvalue at given confidence
	 *
	 * i.e P( | pv - approx | < max_error ) >= confidence
	 *
	 * Beware, the current implementation is very conservative, i.e the actual confidence is much higher
	 * because of the asymmetry of the binomial confidence inervals
	 *
	 * @param the pvalue for which to compute the confidence interval
	 * @param nb_permut the number of permutation used for the approximated pvalue
	 * @param confidence the confidence of the interval
	 *
	 * @return the maximum error
	 */
	static double max_error(double pvalue, INT64 nb_permut,  double confidence);

	/**
	 * compute the maximum relative error for the approximation of the pvalue at a given confidence
	 *
	 * i.e P( | pv - approx | / pv < max_relative_error ) >= confidence
	 *
	 *  see <max_error>
	 *
	 * @param the pvalue for which to compute the confidence interval
	 * @param nb_permut the number of permutation used for the approximated pvalue
	 * @param confidence the confidence of the interval
	 *
	 * @return the maximum error
	 */
	static double max_relative_error(double pvalue, INT64 nb_permut,  double confidence);

	/**
	 * compute the lower bound of an approximated pvalue at a given confidence
	 *
	 * i.e P( pvalue > lower_bound) >= confidence
	 *
	 * @param the pvalue for which to compute the confidence interval
	 * @param nb_permut the number of permutation used for the approximated pvalue
	 * @param confidence the confidence of the interval
	 *
	 * @return true iff the actual pvalue is above
	 */
	static double lower_bound(double pvalue, INT64 nb_permut,  double confidence);

	/**
	 * check if the actual pvalue from an approximated pvalue is above a threshold
	 *
	 * i.e P( pvalue > T) >= confidence
	 *
	 * @param the pvalue for which to compute the confidence interval
	 * @param the threshold
	 * @param nb_permut the number of permutation used for the approximated pvalue
	 * @param confidence the confidence of the interval
	 *
	 * @return true iff the actual pvalue is above
	 */
	static double is_pvalue_above(double pvalue, double threshold, INT64 nb_permut,  double confidence);

private:
	Binomial(); // should prevent instantiation !
};


#endif /* BINOMIAL_HPP_ */
