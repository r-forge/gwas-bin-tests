/*
 * fdr.hpp
 *
 *  Created on: Jul 28, 2010
 *      Author: kforner
 */

#ifndef FDR_HPP_
#define FDR_HPP_

#include "misc.hpp"

class FDR {
private: FDR(); // prevent instantiation

public: // =========== CLASS METHODS ==========
	static void computeRanks(const vector<double>& pvalues, vector<int>& ranks);
	static void computeFDR(const vector<double>& pvalues, vector<double>& fdrs, double pi0 = 1.0);

};

#endif /* FDR_HPP_ */
