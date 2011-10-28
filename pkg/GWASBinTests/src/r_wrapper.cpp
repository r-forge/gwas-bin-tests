#include <Rcpp.h>
#include "nnbc_engine.hpp"
#include "bins.hpp"
#include "chromosomes.hpp"
using namespace Rcpp;

void _split_types(const char* types, vector<string>& types_list) {
	// parse the types string that contain space delimited names

	istringstream in(types);
	string s;
	while ( in >> s)
		types_list.push_back(s);

}


NNBCParameters makeNNBCParams(SEXP params)
{
	// Get parameters in params.
//	RcppParams rparam(params);
	Rcpp::List rparam(params);

	string types_string = as<string>(rparam["types_as_string"]);

	vector<string> types_list;
	_split_types(types_string.c_str(), types_list);
	unsigned long long nb_permutations = as<double>(rparam["nb_permutations"]);
	NNBCParameters parameters(types_list,
			nb_permutations,
			as<int>(rparam["verbosity"]),
			as<int>(rparam["seed"]),
			as<double>(rparam["max_error_rate"]),
			as<int>(rparam["use_affymetrix_model"]),
			as<double>(rparam["min_pvalue"]),
			as<double>(rparam["confidence"]),
			as<double>(rparam["max_relative_error"]),
			as<int>(rparam["regularizeEstimators"]),
			as<int>(rparam["excludeX"]),
			as<int>(rparam["excludeMalesonXChr"])
		);


	return parameters;
}

RcppExport SEXP _runNNBCEngine(NNBCEngine& engine, int nb_threads)
{
	engine.check();

	NNBCResults results;

	// === NNBC ====
//	cerr << "_runNNBCEngine : before process_all_bins\n";
	engine.process_all_bins(results, nb_threads, engine.params().verbose > 0);
	const Bins& thebins = engine.bins();
	const int nb_non_empty_bins = thebins.nb_non_empty_bins();
//	cerr << "nb_non_empty_bins" << nb_non_empty_bins;
	if ( nb_non_empty_bins ) {
		engine.computeFDR(results);

		vector<int> bins(nb_non_empty_bins), chr(nb_non_empty_bins), start(nb_non_empty_bins),
			end(nb_non_empty_bins), nb_snps(nb_non_empty_bins);
		vector<double> nb_permutations_by_bin(nb_non_empty_bins);
		for (int i = 0; i < nb_non_empty_bins; ++i) {
			int bin = thebins.non_empty_bins[i];
			bins[i] = bin + 1; // R starts from 1
			chr[i] = thebins.chrs[bin];
			start[i] = thebins.starts[bin];
			end[i] = thebins.ends[bin];
			nb_snps[i] = thebins.last_snp[bin] - thebins.first_snp[bin] + 1;
		}
		for (int i = 0; i < nb_non_empty_bins; ++i)
			nb_permutations_by_bin[i] = (double)results.nb_permutations_by_bin[i];


		return Rcpp::List::create(
				_["bin"] = bins,
				_["chr"] = chr,
				_["start"] = start,
				_["end"] = end,
				_["nb_snps"] = nb_snps,
				_["nb_permutations"] = nb_permutations_by_bin,
				_["pvalues_by_bin"] = results.pvalues_by_bin,
				_["scores_by_bin"] = results.scores_by_bin,
				_["fdrs_by_bin"] = results.fdrs_by_bin,
				_["types"] = engine.params().getTypes()
		);
	}

	return Rcpp::List::create();



}


RcppExport SEXP RcppRunNNBCOnFiles(SEXP R_files, SEXP params, SEXP groups)
{
	BEGIN_RCPP
	NNBCParameters parameters = makeNNBCParams(params);

	Rcpp::List rparam(params);
	Rcpp::List files(R_files);

	string geno_file = as<string>(files["geno_file"]);
	string pheno_file = as<string>(files["pheno_file"]);
	string bins_file = as<string>(files["bins_file"]);

	vector<int> igroups = as< vector<int> >(groups);

	NNBCEngine engine(geno_file, pheno_file, bins_file, igroups, parameters);

	return _runNNBCEngine(engine, as<int>(rparam["threads"]));

END_RCPP
}


RcppExport SEXP RcppRunNNBCOnTable(SEXP R_nb_phenotypes, SEXP R_vector, SEXP R_params)
{
BEGIN_RCPP

	NNBCParameters parameters = makeNNBCParams(R_params);

	int nb_phenotypes = as<int>(R_nb_phenotypes);
	IntegerVector table_vector(R_vector);
	vector<double> pvalues;
	vector<double> scores;
	int nb_permutations = -1;
	int table_size = SIZE(table_vector) / nb_phenotypes;

	if ( table_size == 4) {  // UnivariateTable
		// make a UnivariateTable
		UnivariateTable T(nb_phenotypes);
		int index = 0;
		// R matrices are column-major
		for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i) {
			for (int p = 0; p < nb_phenotypes; ++p) {
				T(p, i) = table_vector[index++];
			}
		}
		T.fill_for_H0();
		Dataset dataset(T);

		nb_permutations = NNBCEngine::process_dataset(dataset, parameters, pvalues, scores);
	} else if (table_size == 16) {
			// divariate table
			DivariateTable T(nb_phenotypes);
			int index = 0;
			// R matrices are column-major
			for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
				for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j) {
					for (int p = 0; p < nb_phenotypes; ++p) {
						T(p, j, i) = table_vector[index++];
					}
				}
			T.fill_for_H0();
			Dataset dataset(T);
			nb_permutations = NNBCEngine::process_dataset(dataset, parameters, pvalues, scores);
	} else {
		die("Bad table size");
	}

	return Rcpp::List::create(
			_["nb_permutations"] = nb_permutations,
			_["pvalues"] = pvalues,
			_["scores"] = scores,
			_["types"] = parameters.getTypes()
	);

END_RCPP
}


RcppExport SEXP RcppRunNNBConRdata(SEXP idnames
		,SEXP map, SEXP coding, SEXP snpnames, SEXP chromosome, SEXP strand, SEXP gtps_data,
		SEXP phdata, SEXP bins_df, SEXP R_params, SEXP groups	) {
BEGIN_RCPP

	vector<string> chr_names = as< vector<string> >(chromosome);

	// convert it to integer chr codes
	Chromosomes chr_converter;
	vector<int> chrs;
	chrs.reserve( SIZE(chr_names) );
	foreach(string name, chr_names)
		chrs.push_back( chr_converter.convert( name ) );


	GenabelSnpData snp_data(
			as< vector<string> >(idnames)
			,as< vector<int> >(map)
			,as< vector<UBYTE> >(coding)
			,as< vector<string> >(snpnames)
			,chrs
			,as< vector<UBYTE> >(strand)
			,as< vector<UBYTE> >(gtps_data)
		);


	const GenabelSnp& snp0 = snp_data.getSnp(0);

	// PHENOTYPIC DATA CONVERSION - (from data frame)
	DataFrame ph = (phdata);

	// ids : column 0

	vector<string> ids = ph[0];
	vector<UBYTE> sexs = ph[1];
	vector<UBYTE> pops = ph[2];

	vector<int> igroups = as< vector<int> >(groups);
//
//	if ( ph.size() > 3)
//		groups = ph[3];

	GenabelPhenoData pheno_data(ids, sexs, pops);

	// BINS DATA CONVERSION - (from data frame)
	DataFrame bins_tmp(bins_df);
	vector<string> bins_chr_names = bins_tmp[0];
	vector<int> starts = bins_tmp[1];
	vector<int> ends = bins_tmp[2];

//	cerr << "bins_chr_names" << bins_chr_names;
//	cerr << "starts" << starts;
//	cerr << "ends" << ends;

	NNBCParameters parameters = makeNNBCParams(R_params);

	int VERBOSE = parameters.verbose;

	Rcpp::List rparam(R_params);

	Bins bins(bins_chr_names, starts, ends, snp_data, parameters.excludeX);

//	cerr << "bins" << bins;

	NNBCEngine engine(&snp_data, &pheno_data, &bins, igroups, parameters);

	SEXP res= _runNNBCEngine(engine, as<int>(rparam["threads"]));

	return res;

END_RCPP
}

//RcppExport SEXP RcppParamsExample(SEXP params)
//{
//BEGIN_RCPP
//	NNBCParameters parameters = makeNNBCParams(params);
//	cout << parameters;
//END_RCPP
//}

