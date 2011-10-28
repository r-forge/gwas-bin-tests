/*
 * dataset.hpp
 *
 *  Created on: Jun 17, 2010
 *      Author: kforner
 */

#ifndef DATASET_HPP_
#define DATASET_HPP_

#include "misc.hpp"
#include "genabel.hpp"
#include "univariate_table.hpp"
#include "divariate_table.hpp"
//class UnivariateTable;
//class DivariateTable;

class Dataset
{
public:
	Dataset(const vector< vector<UBYTE> >& genotypes, const vector<UBYTE>& phenotypes);
	Dataset(const GenabelSnpData& snp_data, const GenabelPhenoData& pheno_data, int first_snp, int last_snp, bool excludeMales = false);
	// infer one snp
	Dataset(const UnivariateTable& table);
	// infer two snps
	Dataset(const DivariateTable& table);

	/**
	 * make a sub-dataset from an existing dataset by keeping only individuals/samples corresponding to the indices
	 *
	 * @param dataset the dataset to make a sub-dataset from
	 * @param indices the indices of individuals to keep
	 *
	 */
	Dataset(const Dataset& dataset, const vector<int>& indices);

//	Dataset(int nb_snps, int nb_samples, int nb_phenotypes);
//	~Dataset() {}

public: // ========== INSTANCE METHODS ==========
	friend ostream& operator <<(ostream& s, const Dataset& t);

public: // ========== GETTERS ==========
	int nb_phenotypes() const;
	int nb_snps() const;
	int nb_samples() const;

	const vector< vector<UBYTE> >& genotypes() const;
	const vector<UBYTE>& phenotypes() const;
private: // ========== INSTANCE DATA ==========
	int _nb_phenotypes;
	vector< vector<UBYTE> > _genotypes;
	vector<UBYTE> _phenotypes;
};

class DivariateOptimizedDataSet {
public: // ========== LIFECYCLE ==========
	explicit DivariateOptimizedDataSet(Dataset const& dataset);
public: // =========== ACCESSORS ==========
	int nb_phenotypes() const;
	int nb_snps() const;
	int nb_snp_pairs() const;
	int nb_samples() const;
	const vector< vector<UBYTE> >& encodedGenotypePairs() const;
	const vector<UBYTE>& encodedPhenotypes() const;
public: // =========== METHODS ==========
	friend std::ostream& operator<< (std::ostream& o, DivariateOptimizedDataSet const& DivariateOptimizedDataSet);

public: // ========== CLASS METHODS ==========
	static int computePhenotypeShift(int nb_phenotypes);
	static bool isOptimizable(Dataset const& dataset);
	//static void encodeGenotypicPair(vector<UBYTE> const& g1, vector<UBYTE> const& g2, vector<UBYTE>& encoded, int shift);
	static void encodeGenotypicPair(vector<UBYTE> const& g1, vector<UBYTE> const& g2, vector<UBYTE>& encoded);
	static void encodePhenotypes(vector<UBYTE> const& phenotypes, vector<UBYTE>& encoded);

private: // ========== INSTANCE DATA ==========
	int _nb_phenotypes;
	vector< vector<UBYTE> > _encoded_genotypes;
	vector<UBYTE> _encoded_phenotypes;
};

// ========== Dataset IMPLEMENTATION ==========
inline int Dataset::nb_phenotypes() const { return _nb_phenotypes; }
inline int Dataset::nb_snps() const { return _genotypes.size(); }
inline int Dataset::nb_samples() const { return _phenotypes.size(); }
inline const vector< vector<UBYTE> >& Dataset::genotypes() const { return _genotypes; }
inline const vector<UBYTE>& Dataset::phenotypes() const { return _phenotypes; }

// ========== DivariateOptimizedDataSet IMPLEMENTATION ==========
inline int DivariateOptimizedDataSet::nb_phenotypes() const { return _nb_phenotypes; }
inline int DivariateOptimizedDataSet::nb_snp_pairs() const { return _encoded_genotypes.size(); }
inline int DivariateOptimizedDataSet::nb_snps() const { return nb_snp_pairs() + 1; }
inline int DivariateOptimizedDataSet::nb_samples() const { return _encoded_phenotypes.size(); }
inline const vector< vector<UBYTE> >&
DivariateOptimizedDataSet::encodedGenotypePairs() const { return _encoded_genotypes; }

inline const vector<UBYTE>& DivariateOptimizedDataSet::encodedPhenotypes() const { return _encoded_phenotypes; }



#endif /* DATASET_HPP_ */
