/*
 * fdr.cpp
 *
 *  Created on: Jul 28, 2010
 *      Author: kforner
 */

#include "fdr.hpp"



void FDR::computeRanks(const vector<double> & pvalues, vector<int> & ranks)
{
	const int n = SIZE(pvalues);

	// make indices
	vector<int> indices(n, 0);
	for (int i = 0; i < n; ++i)
		indices[i] = i;
	// sort indices

	// WORKAROUND: the following does not work because
	// the var() functor does not accept const arguments
	// std::sort(ALL(indices), var(pvalues)[_1] < var(pvalues)[_2]);
	// so an easy although non-optimal work-around is to duplicate the const pvalues
	// into a non-const one
	vector<double> non_const_pvalues = pvalues;
	std::sort(ALL(indices), var(non_const_pvalues)[_1] < var(non_const_pvalues)[_2]);

	// compute ranks
	ranks.resize(n, 0);
	ranks[indices[n-1]] = n;
	for (int i = n - 2; i >= 0; --i) {
		int index = indices[i];
		int previous = indices[i+1];
		if ( pvalues[index] == pvalues[previous] )
			ranks[index] = ranks[previous];
		else
			ranks[index] = i + 1;
	}
}

void FDR::computeFDR(const vector<double> & pvalues, vector<double> & fdrs, double pi0)
{
	const int n = SIZE(pvalues);

	vector<int> ranks(n);
	computeRanks(pvalues, ranks);

	fdrs.resize(n);
	for (int i = 0; i < n; ++i) {
		fdrs[i] = pvalues[i] * n * pi0 / double(ranks[i]);
	}
}




