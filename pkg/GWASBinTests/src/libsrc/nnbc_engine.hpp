/*
 * nnbc.hpp
 *
 *  Created on: Aug 27, 2010
 *      Author: kforner
 */

#ifndef NNBC_HPP_
#define NNBC_HPP_

#include "misc.hpp"
#include "genabel.hpp"
#include "dataset.hpp"
#include "bins.hpp"

class NNBCParameters {
public:
//	NNBCParameters() {
//		init();
//	}

	/*
	 * Corresponds to the use from command-line nnbc++ executable
	 */

	NNBCParameters(
			const vector<string>& types,
			INT64 nb_permutations=1000,
			int verbosity=1,
			int seed=time(NULL),
			double max_error_rate=0.05,
			bool use_affymetrix_model=false,
			double min_pvalue = -1,
			double confidence = 0.999,
			double max_relative_error = 0,
			bool regularizeEstimators = false,
			bool excludeX = false,
			bool excludeMalesonXChr = false
			);

	void init();
	friend std::ostream& operator<<(std::ostream& o, const NNBCParameters& p);

	vector<string> getTypes() const;

public: // ======= INSTANCE DATA ======
	bool use_affymetrix_model;
	double max_error_rate;
	bool excludeX;
	bool excludeMalesonXChr;
	double regCoeff;
	INT64 nb_permutations;
	int verbose;
	int seed;
	double min_pvalue;
	double confidence;
	double max_relative_error;
	double pi0;
	std::map<string, bool> types;

public: // ===== CONSTANTS =====
	static const string UNIVARIATE;
	static const string DIVARIATE;
	static const string ALLELIC;
	static const string GENOTYPIC;

};

class NNBCResults {
public:
	vector< vector<double> > pvalues_by_bin;
	vector< vector<double> > scores_by_bin;
	vector<INT64> nb_permutations_by_bin;
	vector< vector<double> > fdrs_by_bin;
};

/**
 * Main class to run the NNBC analysis
 */

class NNBCEngine {
public:
	/**
	 * make an engine from a GenABEL dataset
	 *
	 * @param snp_data the data about SNPs (genotypes and genomic locations)
	 * @param pheno_data the phenotypic data
	 * @param bins the information on how to group the SNPs into bins based on their genomic locations
	 *
	 * @params groups an optional vector containing group information for each sample/individual.
	 * 	If supplied it will be used as a covariable in the analysis.
	 * 	If empty and the *gws* column is available in *pheno_data* it will be used as covariable
	 * 	otherwise no covariable is used.
	 *
	 * @param params the NNBC algorithm parameters
	 *
	 *
	 */
	NNBCEngine(GenabelSnpData* snp_data, GenabelPhenoData* pheno_data, 	Bins* bins,
			const vector<int>& groups, const NNBCParameters& params);

	/**
	 * make an engine from a GenABEL dataset stored in files
	 *
	 * @param basename the basename of the files ".gag", ".gap" and ".bins" suffixes will respectively be appended
	 *  	to the basename to form the data file names. Purely for convenience.
	 * @params groups an optional vector containing group information for each sample/individual.
	 * 	If supplied it will be used as a covariable in the analysis.
	 * 	If empty and the *gws* column is available in *pheno_data* it will be used as covariable
	 * 	otherwise no covariable is used.
	 *
	 * @param params the NNBC algorithm parameters
	 *
	 *
	 */
	NNBCEngine(const string& basename, const vector<int>& groups, const NNBCParameters& params);

	/**
	 * make an engine from a GenABEL dataset stored in files
	 *
	 * @param geno_file the GenABEL genotypes file (.gag)
	 * @param pheno_file the GenABEL phenotypes file (.gap)
	 * @param bins_file the bins definition file (.bins)
	 *
	 * @params groups an optional vector containing group information for each sample/individual.
	 * 	If supplied it will be used as a covariable in the analysis.
	 * 	If empty and the *gws* column is available in *pheno_data* it will be used as covariable
	 * 	otherwise no covariable is used.
	 *
	 * @param params the NNBC algorithm parameters
	 *
	 *
	 */
	NNBCEngine(const string&geno_file, const string&pheno_file, const string&bins_file,
			const vector<int>& groups, const NNBCParameters& params);
	~NNBCEngine();

public: // ===== ACCESSORS ======
	const Bins& bins() const { return *_bins; }
	const NNBCParameters& params() const { return _params; }
public:  // ===== METHODS ======
	void check() const;

	INT64 process_bin(int bin, vector<double>& pvalues, vector<double>& scores);

	void process_all_bins(NNBCResults& results, int nb_threads = 0, bool show_progress = false);

	void computeFDR(NNBCResults& results);

	void reportResults(const NNBCResults& results);


public:
	INT64 process_bin(int seed2, const Dataset& dataset, const NNBCParameters& params,
			const GenabelPhenoData& pheno_data, vector<double>& pvalues, vector<double>& scores);

	/**
	 * for group support
	 */
	static INT64 process_bin(int seed2, const Dataset& dataset, const NNBCParameters& params,
			const vector<int>& groups, vector<double>& pvalues, vector<double>& scores);
	/**
	 * no group support
	 */
	static INT64 process_bin(int seed2, const Dataset& dataset, const NNBCParameters& params,
			 vector<double>& pvalues, vector<double>& scores);
//	template <class SHUFFLER>
//		static int process_bin(int seed2, const Dataset& dataset, const NNBCParameters& params
//				, SHUFFLER& shuffler, vector<double>& pvalues, vector<double>& scores);
	static INT64 process_dataset(const Dataset& dataset,
			const NNBCParameters& params, vector<double>& pvalues, vector<double>& scores);

private:
	void init(const vector<int>& groups);
	void load(const string&geno_file, const string&pheno_file, const string&bins_file,  int verbosity = 0, bool excludeXchr = false);


private:  // ======= INSTANCE DATA ======
	NNBCParameters _params;
	GenabelSnpData* _snp_data;
	GenabelPhenoData* _pheno_data;
	Bins* _bins;
	bool _delete_data;
	vector<int> _groups;
};



#endif /* NNBC_HPP_ */
