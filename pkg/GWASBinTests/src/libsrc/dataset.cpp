/*
 * dataset.cpp
 *
 *  Created on: Jun 17, 2010
 *      Author: kforner
 */

#include "dataset.hpp"


ostream& operator <<(ostream& s, const Dataset& d)
{
	s << "Dataset: " << d.nb_snps() << " snps, " << d.nb_samples() << " samples for "
			<< d.nb_phenotypes() << " phenotypes\n";
	return s;
}


Dataset::Dataset(const vector< vector<UBYTE> >& genotypes, const vector<UBYTE>& phenotypes)
	: _genotypes(genotypes), _phenotypes(phenotypes)
{
	_nb_phenotypes = *max_element(phenotypes.begin(), phenotypes.end()) + 1;
}

Dataset::Dataset(const GenabelSnpData& snp_data,
		const GenabelPhenoData& pheno_data, int first_snp, int last_snp, bool excludeMales)
	:  _nb_phenotypes(-1), _genotypes(), _phenotypes(pheno_data.getPops())
{
	const int nb_snps = last_snp - first_snp + 1;
	const int nb_samples = snp_data.nb_samples();
	const int MALE = GenabelPhenoData::MALE; // instantiate the constant
	const vector<UBYTE> sexs = pheno_data.getSexs();
	int nb_samples_to_keep = nb_samples;
	if ( excludeMales ) {
		// count males
		int males = std::count( ALL(sexs), MALE );
		nb_samples_to_keep -= males;
	}
	_genotypes.resize(nb_snps);
	vector<UBYTE> genotypes(nb_samples);
	for (int i = 0; i < nb_snps; ++i) {
		int snp_idx = i + first_snp;

		snp_data.getSnp(snp_idx).decodeGenotypes( genotypes );
		if ( excludeMales ) {
			_genotypes[i].reserve(nb_samples_to_keep);
			for (int j = 0; j < nb_samples; ++j) {
				if ( sexs[j] != MALE )
					_genotypes[i].push_back( genotypes[j] );
			}
			assert( SIZE(_genotypes[i]) == nb_samples_to_keep );
		} else {
			_genotypes[i] = genotypes;
		}
	}
	_nb_phenotypes = *max_element(_phenotypes.begin(), _phenotypes.end()) + 1;
}


Dataset::Dataset(const UnivariateTable& table)
	: _genotypes(1), _phenotypes()
{
	const int nb_samples = table.total(true);
	const int nb_phenotypes = table.nb_phenotypes();
	_nb_phenotypes = nb_phenotypes;

	vector<UBYTE>& genos = _genotypes[0];
	genos.resize(nb_samples, 0);
	_phenotypes.resize(nb_samples, 0);
	int cursor = 0;
	for (int k = 0; k < nb_phenotypes; ++k) {
		for (int g = 0; g < NB_GENOTYPES_WITH_ZZ; ++g) {
			const int nb = table[k][g];
			for (int i = 0; i < nb; ++i) {
				genos[cursor] = g;
				_phenotypes[cursor] = k;
				cursor++;
			}
		}
	}

	assert( cursor == nb_samples );
}

Dataset::Dataset(const DivariateTable& table)
	: _genotypes(2), _phenotypes()
{
	const int nb_samples = table.total(true);
	const int nb_phenotypes = table.nb_phenotypes();
	_nb_phenotypes = nb_phenotypes;

	vector<UBYTE>& g0 = _genotypes[0];
	vector<UBYTE>& g1 = _genotypes[1];
	g0.resize(nb_samples, 0);
	g1.resize(nb_samples, 0);
	_phenotypes.resize(nb_samples, 0);
	int cursor = 0;
	for (int k = 0; k < nb_phenotypes; ++k) {
		for (int g = 0; g < NB_GENOTYPES_WITH_ZZ; ++g)
			for (int g2 = 0; g2 < NB_GENOTYPES_WITH_ZZ; ++g2) {
				const int nb = table[k][g][g2];
				for (int i = 0; i < nb; ++i) {
					g0[cursor] = g;
					g1[cursor] = g2;
					_phenotypes[cursor] = k;
					cursor++;
				}
			}
	}

	assert( cursor == nb_samples );
}

Dataset::Dataset(const Dataset& dataset, const vector<int>& indices)
	: _nb_phenotypes(dataset.nb_phenotypes()), _genotypes(dataset.nb_snps()), _phenotypes(SIZE(indices))
{
	int new_size = SIZE(indices);
	int nb_snps = dataset.nb_snps();

	const vector<UBYTE>& source_phenotypes = dataset.phenotypes();
	const vector< vector<UBYTE> >& source_genotypes = dataset.genotypes();

	for (int i = 0; i < new_size; ++i)
		_phenotypes[i] = source_phenotypes[indices[i]];

	for (int snp = 0; snp < nb_snps; ++snp) {
		const vector<UBYTE>& source_genos = source_genotypes[snp];
		vector<UBYTE>& genos = _genotypes[snp];
		genos.resize(new_size);
		for (int i = 0; i < new_size; ++i)
			genos[i] = source_genos[indices[i]];
	}
}


