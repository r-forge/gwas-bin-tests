/*
 * binomial.cpp
 *
 *  Created on: Sep 20, 2011
 */
#include "binomial.hpp"



Interval Binomial::confidence_interval(double pvalue, INT64 nb_permut, double confidence)
{
	const double alpha = (1.0 - confidence) / 2.0;

	double success = pvalue * nb_permut;
	Interval i;

	i.first = binomial_distribution<>::find_lower_bound_on_p(nb_permut, success, alpha);
	i.second = binomial_distribution<>::find_upper_bound_on_p(nb_permut, success, alpha);

	return i;
}

double Binomial::max_error(double pvalue, INT64 nb_permut, double confidence)
{
	const double alpha = (1.0 - confidence) / 2.0;
	double success = pvalue * nb_permut;
	double lower = binomial_distribution<>::find_lower_bound_on_p(nb_permut, success, alpha);
	double upper = binomial_distribution<>::find_upper_bound_on_p(nb_permut, success, alpha);
	return max(pvalue-lower, upper-pvalue);

}

double Binomial::max_relative_error(double pvalue, INT64 nb_permut, double confidence)
{
	return max_error(pvalue, nb_permut, confidence) / pvalue;
}

double Binomial::lower_bound(double pvalue, INT64 nb_permut,  double confidence)
{
	return binomial_distribution<>::find_lower_bound_on_p(nb_permut, pvalue*nb_permut, 1-confidence);
}

double Binomial::is_pvalue_above(double pvalue, double threshold, INT64 nb_permut,  double confidence)
{
	return lower_bound(pvalue, nb_permut, confidence) >= threshold;
}
