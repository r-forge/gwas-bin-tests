/*
 * divariateTable.cpp
 *
 *  Created on: Jun 23, 2010
 *      Author: kforner
 */

#include "divariate_table.hpp"
#include <cmath>
#include "dataset.hpp"

DivariateTable::DivariateTable(int nb_phenotypes) :
		_nb_true_phenotypes(nb_phenotypes),
		_nb_phenotypes_with_h0(nb_phenotypes+1),
		_h0_phenotype(nb_phenotypes),
		_matrix(nb_phenotypes+1, NB_GENOTYPES_WITH_ZZ, NB_GENOTYPES_WITH_ZZ)
{}

bool DivariateTable::operator ==(DivariateTable const& t) const {
	return _matrix == t._matrix;
}

void DivariateTable::fill(const UBYTE* genotypes1, const UBYTE* genotypes2
		, const UBYTE *phenotypes, int size)
{
	Matrix3D<int>& m = _matrix;
	m.fill(0); // important

	for (int row = 0; row < size; ++row) {
		m(phenotypes[row], genotypes1[row], genotypes2[row])++;
	}

	fill_for_H0();
}



void DivariateTable::fill_for_H0() {
	const int nb_phenotypes = _nb_true_phenotypes;
	int** h0 = _matrix[_h0_phenotype];
	// MAYBE: try to replace by a single loop
	for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
		for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
			h0[i][j] = 0;
	int** row;
	for (int k = 0; k < nb_phenotypes; ++k) {
		row = _matrix[k];
		// MAYBE: try to replace by a single loop
		for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
			for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
				h0[i][j] += row[i][j];
	}
}

void DivariateTable::fill_nphenotypes(const UBYTE *encoded_genotype_pairs, const UBYTE *encoded_phenotypes, int size)
{
	if ( nb_phenotypes() == 2) {
		fill_2phenotypes(encoded_genotype_pairs, encoded_phenotypes, size);
		return;
	}

	// the way the genotypes were encoded should exactly match
	// the Matrix vector order

	_matrix.fill(0);
	int* T = _matrix.c_contiguous_array();

	for (int i = 0; i < size; ++i)
		T[encoded_genotype_pairs[i] | encoded_phenotypes[i]]++;

	fill_for_H0();
}

void DivariateTable::fill_2phenotypes(const UBYTE *encoded_genotype_pairs,
		const UBYTE *encoded_phenotypes, int size)
{
	int*** m = _matrix.x_begin();
	int* T = &m[0][0][0];
	int** t0 = m[0];
	int** t1 = m[1];
	int** h0 = m[_h0_phenotype];

    for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
		for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
			t0[i][j] = t1[i][j] = 0;

    // all the computing time is spent in this loop
	for (int row = 0; row < size; ++row)
		T[encoded_genotype_pairs[row] | encoded_phenotypes[row]]++;

    for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
		for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
			h0[i][j] = t0[i][j] + t1[i][j];

}


std::ostream& operator<< (std::ostream& o,  const DivariateTable& t) {
	return o << t._matrix;
}




DivariateOptimizedDataSet::DivariateOptimizedDataSet(const Dataset & dataset) :
	_nb_phenotypes(dataset.nb_phenotypes()), _encoded_genotypes(
			dataset.nb_snps() - 1), _encoded_phenotypes(dataset.phenotypes()) {

	// check it: nb_phenotypes must be codable on 4 bits ie <= 16
	if ( ! isOptimizable(dataset) ) {
		cerr << "Too many phenotypes for DivariateOptimizedDataSet !! " << dataset.nb_phenotypes() << "\n";
		die("Too many phenotypes for DivariateOptimizedDataSet !!" );
	}
	assert( isOptimizable(dataset) );

	encodePhenotypes(dataset.phenotypes(), _encoded_phenotypes);
//	const int shift = computePhenotypeShift(nb_phenotypes());

	const int nb = nb_snp_pairs();
	for (int i = 0; i < nb; ++i) {
//		encodeGenotypicPair( dataset.genotypes()[i], dataset.genotypes()[i+1],
//				_encoded_genotypes[i], shift);
		encodeGenotypicPair( dataset.genotypes()[i], dataset.genotypes()[i+1],
				_encoded_genotypes[i]);
	}
}

int DivariateOptimizedDataSet::computePhenotypeShift(int nb_phenotypes)
{
	return (int)ceil( log2(nb_phenotypes) );
}

void DivariateOptimizedDataSet::encodeGenotypicPair(
		vector<UBYTE> const& g1, vector<UBYTE> const& g2,
		vector<UBYTE>& encoded) {
	const int nb = SIZE(g1);
	assert( SIZE(g2) == nb);
	encoded.resize(nb, 0);
	for (int i = 0; i < nb; ++i) {
		encoded[i] = (g2[i] + (g1[i] << 2) );
	}
}

void DivariateOptimizedDataSet::encodePhenotypes(
		vector<UBYTE> const& phenotypes, vector<UBYTE>& encoded) {
	const int nb = SIZE(phenotypes);
	encoded.resize(nb, 0);
	for (int i = 0; i < nb; ++i) {
		encoded[i] = phenotypes[i] << 4;
	}
}



//void DivariateOptimizedDataSet::encodeGenotypicPair(
//		vector<UBYTE> const& g1, vector<UBYTE> const& g2,
//		vector<UBYTE>& encoded, int shift) {
//	const int nb = SIZE(g1);
//	assert( SIZE(g2) == nb);
//	encoded.resize(nb, 0);
//	for (int i = 0; i < nb; ++i) {
//		encoded[i] = (g2[i] + (g1[i] << 2) ) << shift;
//	}
//}

bool DivariateOptimizedDataSet::isOptimizable(const Dataset & dataset)
{
	return dataset.nb_phenotypes()  <= (1<<4);
}

std::ostream & operator <<(std::ostream & o, const DivariateOptimizedDataSet & d)
{
	o << "DivariateOptimizedDataSet: " << d.nb_snps() << " snps ("<< d.nb_snp_pairs() << " pairs), "
			<< d.nb_samples() << " samples for "
			<< d.nb_phenotypes() << " phenotypes\n";
	return o;
}


