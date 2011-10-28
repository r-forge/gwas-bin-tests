/*
 * pvalues_computable.hpp
 *
 *  Created on: Sep 20, 2011
 *
 */



#ifndef PVALUES_COMPUTABLE_HPP_
#define PVALUES_COMPUTABLE_HPP_

#include "misc.hpp"
#include "shuffle.hpp"

 /**
 * This class is an abstract class for implementing pvalue computation for score statistics.
 * It is meant to be use by a <PvaluesComputer> instance.
 *
 * To implement your own p-value computable score statistics, just subclass PvaluesComputable
 * and implement and override if needed the main interface methods;
 *
 */
class PvaluesComputable {

public: // ==== MAIN INTERFACE TO IMPLEMENT =====
	//virtual void init() = 0;
	virtual void compute(vector<double>& scores) = 0;
	virtual void shuffle() = 0;

public: // ===== LIFECYCLE =====
	virtual ~PvaluesComputable() {}

};

/**
* This class is an helper class for statistics scores whose pvalues
* are based on phenotypes shuffling. It currently uses the StdShuffler.
*
* All you have to do is to subclass it and implement the compute(phenotypes, scores) method
*
*/
class PhenotypesShufflingPvaluesComputable : public PvaluesComputable {

public: // ==== PhenotypesShufflingPvaluesComputable INTERFACE TO IMPLEMENT =====
	virtual void compute_for_phenotypes(const vector<UBYTE>& phenotypes, vector<double>& scores) = 0;
public: // ===== LIFECYCLE =====
	PhenotypesShufflingPvaluesComputable(const vector<UBYTE>& phenotypes, StdShuffler& shuffler)
			: _phenotypes(phenotypes), _shuffler(shuffler) {}
//	PhenotypesShufflingPvaluesComputable(const vector<UBYTE>& phenotypes, int seed1, int seed2)
//		: _phenotypes(phenotypes), _shuffler(StdShuffler(seed1, seed2)) {}
	virtual ~PhenotypesShufflingPvaluesComputable() {}

public: // ==== IMPLEMENTATION OF PvaluesComputable MAIN INTERFACE =====
	//virtual void init() {}
	virtual void shuffle() { _shuffler.shuffle(_phenotypes); }
	virtual void compute(vector<double>& scores) { compute_for_phenotypes(_phenotypes, scores); }

public: // ==== ACCESSORS/SETTERS =====
	void set_phenotypes(const vector<UBYTE>& phenotypes) { _phenotypes = phenotypes; }
	const vector<UBYTE>& get_phenotypes() const { return _phenotypes; }

private: // ==== ATTRIBUTES ====
	vector<UBYTE> _phenotypes;
	StdShuffler& _shuffler;
};


#endif /* PVALUES_COMPUTABLE_HPP_ */
