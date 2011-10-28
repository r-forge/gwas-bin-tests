

/*
 * nnbc.cpp
 *
 *  Created on: Aug 27, 2010
 *      Author: kforner
 */

#ifdef _OPENMP
#include <omp.h>
#endif

//#ifdef USING_R
//void R_CheckUserInterrupt(void);
//#endif

#include <boost/progress.hpp>
#include <boost/timer.hpp>

#include "nnbc_engine.hpp"
#include "genotype_error_model.hpp"
//#include "pvalues.hpp"
#include "iterative_table_scores_pvalues_computable.hpp"
#include "group_iterative_table_scores_pvalues_computable.hpp"
#include "univariate_score_function.hpp"
#include "divariate_score_function.hpp"
#include "sum_score_functions.hpp"
#include "fdr.hpp"
#include "chromosomes.hpp"
#include "pvalues_computer.hpp"

const string NNBCParameters::UNIVARIATE = "univariate";
const string NNBCParameters::DIVARIATE = "divariate";
const string NNBCParameters::ALLELIC = "allelic";
const string NNBCParameters::GENOTYPIC = "genotypic";

NNBCParameters::NNBCParameters(const vector<string>& types_list,
		INT64 nb_permutations, int verbosity,
		int seed, double max_error_rate,
		bool use_affymetrix_model, double min_pvalue,
		double confidence, double max_relative_error,
		bool regularizeEstimators,
		bool excludeX,
		bool excludeMalesonXChr)
{
	init();

	this->nb_permutations = nb_permutations;
	this->verbose = verbosity;
	this->seed = seed;
	this->max_error_rate = max_error_rate;
	this->use_affymetrix_model = use_affymetrix_model;
	this->min_pvalue = min_pvalue;
	this->confidence = confidence;
	this->max_relative_error = max_relative_error;

	if ( regularizeEstimators )
		this->regCoeff = this->max_error_rate;

	vector<string> TYPES;
	TYPES += UNIVARIATE, DIVARIATE, ALLELIC, GENOTYPIC;
	std::map<string, bool> types_map;
	foreach(string t, TYPES) { types_map[t] = true; }

	if (types_list.empty())
		foreach(string t, TYPES)
			{ this->types[t] = true;}
	else
		foreach(string t, types_list) {
			if ( ! isInMap(types_map, t) ) {
				cerr << "Bad type: [" << t << "]\n";
				cerr << "type must be in " << TYPES;
				die("Bad Parameter");
			}
			this->types[t] = true;

		}
}

vector<string> NNBCParameters::getTypes() const
{
	vector<string> TYPES;
	TYPES += UNIVARIATE, DIVARIATE, ALLELIC, GENOTYPIC;
	vector<string> res;
	foreach(string t, TYPES)
		if ( isInMap(types, t) )
			res.push_back(t);
	return res;
}

void NNBCParameters::init()
{
	use_affymetrix_model = false;
	max_error_rate = 0.05;
	excludeX = false;
	excludeMalesonXChr = false;
	regCoeff = 0;
	nb_permutations = 1000;
	verbose = 0;
	seed = 0;
	min_pvalue = -1;
	confidence = 0.999;
	max_relative_error = -1;
	pi0 = 1.0;
}

std::ostream& operator<<(std::ostream& o, const NNBCParameters& p)
{
	o << "NNBC Parameters:\n";
	o << "use_affymetrix_model=" << p.use_affymetrix_model << "\n";
	o << "excludeX=" << p.excludeX << "\n";
	o << "excludeMalesonXChr=" << p.excludeMalesonXChr << "\n";
	o << "regCoeff=" << p.regCoeff << "\n";
	o << "nb_permutations=" << p.nb_permutations << "\n";
	o << "verbose=" << p.verbose << "\n";
	o << "seed=" << p.seed << "\n";
	o << "min_pvalue=" << p.min_pvalue << "\n";
	o << "confidence=" << p.confidence << "\n";
	o << "max_relative_error=" << p.max_relative_error << "\n";
	o << "pi0=" << p.pi0 << "\n";
	o << "types=\n";
	vector<string> TYPES;
	TYPES += NNBCParameters::UNIVARIATE, NNBCParameters::DIVARIATE,
			NNBCParameters::ALLELIC, NNBCParameters::GENOTYPIC;

	foreach(string s, TYPES)
		if (isInMap(p.types, s))
			o << "\t" << s << "\n";

	return o;
}

void NNBCEngine::load(const string & geno_file, const string & pheno_file, const string & bins_file, int verbosity, bool excludeXchr)
{
	_snp_data = 0;
	_pheno_data = 0;
	_bins = 0;
	_delete_data = true;

	// ==== load data ====
	_snp_data = new GenabelSnpData(geno_file);
	if (verbosity > 0)
	cerr << "loaded snp data: " << "nb_snps=" << _snp_data->nb_snps()
	<< ", nb_samples=" << _snp_data->nb_samples() << endl;

	_pheno_data = new GenabelPhenoData(pheno_file);
	if (verbosity > 0)
	cerr << "loaded pheno_data\n";
	_bins = new Bins(bins_file.c_str(), *_snp_data, excludeXchr );
	if (verbosity > 0)
	cerr << "loaded bins: " << *_bins;

	// ========================================


}

NNBCEngine::NNBCEngine(GenabelSnpData* snp_data, GenabelPhenoData* pheno_data, 	Bins* bins,
		const vector<int>& groups, const NNBCParameters& params)
	: _params(params), _snp_data(snp_data), _pheno_data(pheno_data), _bins(bins), _delete_data(false)
{
	init(groups);
}

NNBCEngine::NNBCEngine(const string & basename, const vector<int>& groups, const NNBCParameters& params)
 : _params(params), _delete_data(true)

{
	string geno_file = basename + ".raw";
	string pheno_file = basename + ".dat";
	string bins_file = basename + ".bins";
	load(geno_file, pheno_file, bins_file, _params.verbose);
	init(groups);
}

NNBCEngine::NNBCEngine(const string & geno_file, const string & pheno_file, const string & bins_file,
		const vector<int>& groups, const NNBCParameters& params)
: _params(params), _delete_data(true)
{
	load(geno_file, pheno_file, bins_file, _params.verbose);
	init(groups);
}

void NNBCEngine::init(const vector<int>& groups)
{
	if ( ! groups.empty() ) {
		_groups = groups;
	} else {
		if ( _pheno_data->max_groups() > 1) {
			if ( _params.verbose > 0)
				REprintf("Using gws(study) variable as default co-variable \n");
			_groups = _pheno_data->getGroups();
		}
	}
	check();
}

NNBCEngine::~NNBCEngine()
{
	if ( _delete_data ) {
		delete _snp_data;
		delete _pheno_data;
		delete _bins;
	}
}

void NNBCEngine::process_all_bins(NNBCResults& results, int nb_threads, bool show_progress)
{
#ifdef _OPENMP
	if ( nb_threads > 0 )
	omp_set_num_threads( nb_threads );
	if ( params().verbose )
		cerr << "Running under OPENMP with " << omp_get_max_threads() << " thread(s)\n";
#endif

	const int nb_non_empty_bins = bins().nb_non_empty_bins();

	if ( params().verbose ) {
		cerr << "nb_non_empty_bins=" << nb_non_empty_bins << endl;

	}

	results.pvalues_by_bin.resize(nb_non_empty_bins);
	results.scores_by_bin.resize(nb_non_empty_bins);
	results.nb_permutations_by_bin.resize(nb_non_empty_bins);

	boost::progress_display* progress = 0;
	if ( show_progress )
		progress = new boost::progress_display( nb_non_empty_bins, cerr );


	// stopping a OMP parallel loop, cf http://www.thinkingparallel.com/2007/06/29/breaking-out-of-loops-in-openmp/
	//bool abort = false;

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < nb_non_empty_bins; ++i) {
		if ( show_progress )
			++*progress;

		int bin = bins().non_empty_bins[i];
		results.nb_permutations_by_bin[i] = process_bin(bin, results.pvalues_by_bin[i], results.scores_by_bin[i]);
	}

	if ( show_progress != 0) {
		cerr << endl;
		delete progress;
	}
}

void NNBCEngine::computeFDR(NNBCResults& results)
{
	const int nb_non_empty_bins = bins().nb_non_empty_bins();
	vector<vector<double> >& fdrs_by_bin = results.fdrs_by_bin;
	const vector< vector<double> >& pvalues_by_bin = results.pvalues_by_bin;
	const int nb_pvalues = pvalues_by_bin.front().size();

	fdrs_by_bin.resize(nb_non_empty_bins);
	for (int i = 0; i < nb_non_empty_bins; ++i)
		fdrs_by_bin[i].resize(nb_pvalues);

	vector<double> pvalues(nb_non_empty_bins);
	vector<double> fdrs(nb_non_empty_bins);
	for (int p = 0; p < nb_pvalues; ++p)
	{
		for (int i = 0; i < nb_non_empty_bins; ++i)
			pvalues[i] = pvalues_by_bin[i][p];
		FDR::computeFDR(pvalues, fdrs, _params.pi0);
		for (int i = 0; i < nb_non_empty_bins; ++i)
			fdrs_by_bin[i][p] = fdrs[i];
	}
}

void NNBCEngine::reportResults(const NNBCResults& results)
{
	const vector<vector<double> >& fdrs_by_bin = results.fdrs_by_bin;
	const vector< vector<double> >& pvalues_by_bin = results.pvalues_by_bin;

	const int nb_pvalues = pvalues_by_bin.front().size();

	const char delim = '\t';
	string PV_UNIV = "pv_univ";
	string PV_DIV = "pv_div";
	vector<string> headers;
	headers += "bin", "chr", "start", "end", "nb_snps", "nb_permuts";

	vector<string> TYPES;
	TYPES += NNBCParameters::UNIVARIATE, NNBCParameters::DIVARIATE,
			NNBCParameters::ALLELIC, NNBCParameters::GENOTYPIC;

	foreach(string s, TYPES)
	if (_params.types[s]) {
		headers += ("pv_" + s);
		headers += ("score_" + s);
		if ( ! fdrs_by_bin.empty() )
			headers += ("FDR_" + s);
	}

	cout << headers.front();
	for_each(headers.begin()+1, headers.end(), cout << constant(delim) << _1);
	cout << endl;
	vector<double> fdrs;
	const Bins& thebins = bins();

	const int nb_non_empty_bins = thebins.nb_non_empty_bins();
	for (int i = 0; i < nb_non_empty_bins; ++i) {
		int bin = thebins.non_empty_bins[i];
		cout << (bin+1) << delim << thebins.chrs[bin] << delim << thebins.starts[bin] << delim
		<< thebins.ends[bin] << delim << (thebins.last_snp[bin] - thebins.first_snp[bin] + 1)
		<< delim << results.nb_permutations_by_bin[i];

		for (int p = 0; p < nb_pvalues; ++p) {
			cout << delim << pvalues_by_bin[i][p];
			cout << delim << results.scores_by_bin[i][p];
			if ( ! fdrs_by_bin.empty() )
				cout << delim << fdrs_by_bin[i][p];
		}

		cout << "\n";
	}
}

void NNBCEngine::check() const
{
	_snp_data->check();
	_pheno_data->check();
	_bins->check();
	if ( _snp_data->nb_samples() != _pheno_data->nb_samples() )
		REprintf("check failed in NNBCEngine: bad nb_samples between snp_data and pheno_data\n");

	if ( !_groups.empty() && SIZE(_groups) != _pheno_data->nb_samples()) {
		REprintf("check failed in NNBCEngine: bad groups length \n");
	}
}

INT64 NNBCEngine::process_bin(int bin, vector<double>& pvalues, vector<double>& scores)
{
	int start = bins().first_snp[bin];
	int end = bins().last_snp[bin];
	Dataset dataset(*_snp_data, *_pheno_data, start, end,
			_params.excludeMalesonXChr && bins().chrs[bin] == Chromosomes::X);

	return ( ! _groups.empty() ) ?
			process_bin(bin, dataset, _params, _groups, pvalues, scores)
			:
			process_bin(bin, dataset, _params, pvalues, scores);

}

INT64 NNBCEngine::process_bin(int seed2, const Dataset & dataset, const NNBCParameters & params,
		 vector<double> & pvalues, vector<double>& scores)
{
	vector<int> groups;
	return process_bin(seed2, dataset, params, groups, pvalues, scores);
}

INT64 NNBCEngine::process_bin(int seed2, const Dataset & dataset, const NNBCParameters & params,
		const vector<int>& groups, vector<double> & pvalues, vector<double>& scores)
{
	GenotypeErrorModel* error_model = 0;
	if (params.use_affymetrix_model)
		error_model = new AffymetrixGenotypeErrorModel(dataset.genotypes(),
				params.max_error_rate);
	else {
		error_model = new SimpleGenotypeErrorModel(dataset.genotypes(),
				params.max_error_rate);
		if ( params.verbose > 1) {
			cerr << "SimpleGenotypeErrorModel: alpha=" << params.max_error_rate << ", beta=" << error_model->getErrorModelForSnp(0)[3][0] << endl;
		}
	}

	const int nb_phenotypes = dataset.nb_phenotypes();

	vector<IterativeTableScoreFunction*> score_functions;
	if (isInMap(params.types, NNBCParameters::UNIVARIATE))
		score_functions.push_back(new UnivariateScoreFunction(nb_phenotypes,
				params.regCoeff, *error_model));
	if (isInMap(params.types, NNBCParameters::DIVARIATE))
		score_functions.push_back(new DivariateScoreFunction(nb_phenotypes,
				params.regCoeff, *error_model));
	if (isInMap(params.types, NNBCParameters::ALLELIC))
		score_functions.push_back(new AllelicSumScoreFunction(nb_phenotypes,
				dataset.nb_snps()));
	if (isInMap(params.types, NNBCParameters::GENOTYPIC))
		score_functions.push_back(new GenotypicSumScoreFunction(nb_phenotypes,
				dataset.nb_snps()));

	INT64 nb_permutations = 0;
	PvaluesComputer pv_computer;
	pv_computer.set_confidence(params.confidence);
	pv_computer.set_min_pvalue(params.min_pvalue);
	pv_computer.set_max_relative_error(params.max_relative_error);

	StdShuffler shuffler(params.seed, seed2);

	PvaluesComputable* comp = 0;
	if ( ! groups.empty() ) {
		comp = new GroupIterativeTableScoresPvaluesComputable(score_functions, dataset, groups, shuffler);
	} else
		 comp = new IterativeTableScoresPvaluesComputable(score_functions, dataset, shuffler);

	nb_permutations = pv_computer.computePvalues(params.nb_permutations, *comp);

	scores = pv_computer.get_observed_scores();
	pvalues = pv_computer.get_pvalues();

	// cleanup
	delete error_model;
	const int nb = SIZE(score_functions);
	for (int i = 0; i < nb; ++i)
		delete score_functions[i];
	delete comp;

	return nb_permutations;



//	if (pheno_data.max_groups() > 1)
//	{
//		const vector<UBYTE>& g = pheno_data.getGroups();
//		StdGroupShuffler shuffler(g, params.seed, seed2);
//		return process_bin(seed2, dataset, params, shuffler, pvalues, scores);
//	}
//	else
//	{
//		StdShuffler shuffler(params.seed, seed2);
//		return process_bin(seed2, dataset, params, shuffler, pvalues, scores);
//	}
}

//template <class SHUFFLER>
//int NNBCEngine::process_bin(int seed2, const Dataset& dataset, const NNBCParameters& params
//		, SHUFFLER& shuffler, vector<double>& pvalues, vector<double>& scores)
//{
//
//	GenotypeErrorModel* error_model = 0;
//	if (params.use_affymetrix_model)
//		error_model = new AffymetrixGenotypeErrorModel(dataset.genotypes(),
//				params.max_error_rate);
//	else {
//		error_model = new SimpleGenotypeErrorModel(dataset.genotypes(),
//				params.max_error_rate);
//		if ( params.verbose > 1) {
//			cerr << "SimpleGenotypeErrorModel: alpha=" << params.max_error_rate << ", beta=" << error_model->getErrorModelForSnp(0)[3][0] << endl;
//		}
//	}
//
//	const int nb_phenotypes = dataset.nb_phenotypes();
//
//	vector<IterativeTableScoreFunction*> score_functions;
//	if (isInMap(params.types, NNBCParameters::UNIVARIATE))
//		score_functions.push_back(new UnivariateScoreFunction(nb_phenotypes,
//				params.regCoeff, *error_model));
//	if (isInMap(params.types, NNBCParameters::DIVARIATE))
//		score_functions.push_back(new DivariateScoreFunction(nb_phenotypes,
//				params.regCoeff, *error_model));
//	if (isInMap(params.types, NNBCParameters::ALLELIC))
//		score_functions.push_back(new AllelicSumScoreFunction(nb_phenotypes,
//				dataset.nb_snps()));
//	if (isInMap(params.types, NNBCParameters::GENOTYPIC))
//		score_functions.push_back(new GenotypicSumScoreFunction(nb_phenotypes,
//				dataset.nb_snps()));
//
//
//	int nb_permutations = 0;
//	PvaluesComputer pv_computer;
//	pv_computer.set_confidence(params.confidence);
//	pv_computer.set_min_pvalue(params.min_pvalue);
//	pv_computer.set_max_relative_error(params.max_relative_error);
//
//
////	IterativeTableScoresPvaluesComputable computable(score_functions, dataset, params.seed, seed2);
//	StdShuffler shuf( params.seed, seed2);
//	IterativeTableScoresPvaluesComputable computable(score_functions, dataset, shuf);
//
//	nb_permutations = pv_computer.computePvalues(params.nb_permutations, computable);
//
//
//	delete error_model;
//	const int nb = SIZE(score_functions);
//	for (int i = 0; i < nb; ++i)
//		delete score_functions[i];
//
//	scores = pv_computer.get_observed_scores();
//	pvalues = pv_computer.get_pvalues();
//	return nb_permutations;
//}

INT64 NNBCEngine::process_dataset(const Dataset& dataset,
			const NNBCParameters& params, vector<double>& pvalues, vector<double>& scores)
{
	StdShuffler shuffler(params.seed);
	return process_bin(0, dataset, params, pvalues, scores);
}

