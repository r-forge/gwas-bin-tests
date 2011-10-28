

/*
 * divariate_table.cpp
 *
 *  Created on: Jun 22, 2010
 *      Author: kforner
 */

#include "table2snp.hpp"
#include "dataset.hpp"

Table2SNP::Table2SNP(int nb_phenotypes)
{
	assert(nb_phenotypes >= 2);
	if ( nb_phenotypes == 2)
		_table = new Table2SNPImpl2Phenotypes();
	else
		_table = new Table2SNPImplGeneric(nb_phenotypes);
}



std::ostream& operator<< (std::ostream& o,  const Table2SNP& t) {
	return t._table->printOn(o);
}

Table2SNP::~Table2SNP()
{
	delete _table;
}

Table2SNPImplGeneric::Table2SNPImplGeneric(int nb_phenotypes) :
		_nb_true_phenotypes(nb_phenotypes),
		_nb_phenotypes_with_h0(nb_phenotypes+1),
		_h0_phenotype(nb_phenotypes),
		_matrix(nb_phenotypes+1, NB_GENOTYPES_WITH_ZZ, NB_GENOTYPES_WITH_ZZ)
{}

Table2SNPImplGeneric::~Table2SNPImplGeneric() {}

void Table2SNPImplGeneric::fill(const UBYTE* genotypes1, const UBYTE* genotypes2
		, const UBYTE *phenotypes, int size)
{
	Matrix3D<int>& m = _matrix;
	m.fill(0); // important

	for (int row = 0; row < size; ++row) {
		m(phenotypes[row], genotypes1[row], genotypes2[row])++;
	}

	fill_for_H0();
}

void Table2SNPImplGeneric::fill_for_H0() {
	const int nb_phenotypes = _nb_true_phenotypes;
	int** h0 = _matrix[_h0_phenotype];
	for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
		for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
			h0[i][j] = 0;
	int** row;
	for (int k = 0; k < nb_phenotypes; ++k) {
		row = _matrix[k];
		for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
			for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
				h0[i][j] += row[i][j];
	}
}

void Table2SNPImplGeneric::fill(const UBYTE *encoded_genotype_pairs, const UBYTE *encoded_phenotypes, int size)
{
	// the way the genotypes were encoded should exactly match
	// the Matrix vector order

	_matrix.fill(0);
	int* T = _matrix.c_contiguous_array();

	for (int i = 0; i < size; ++i)
		T[encoded_genotype_pairs[i] | encoded_phenotypes[i]]++;

	fill_for_H0();
}

ostream & Table2SNPImplGeneric::printOn(std::ostream & o) const
{
	return o << _matrix;
}


void Table2SNPImpl2Phenotypes::fill(const UBYTE* genotypes1, const UBYTE* genotypes2
		, const UBYTE *phenotypes, int size)
{
	int (*m)[NB_GENOTYPES_WITH_ZZ][NB_GENOTYPES_WITH_ZZ] = _matrix;
	int (*t0)[NB_GENOTYPES_WITH_ZZ] = m[0];
	int (*t1)[NB_GENOTYPES_WITH_ZZ] = m[1];
	int (*h0)[NB_GENOTYPES_WITH_ZZ] = m[_H0_PHENOTYPE];


    for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
		for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
			t0[i][j] = t1[i][j] = 0;

	for (int row = 0; row < size; ++row)
		m[phenotypes[row]][genotypes1[row]][genotypes2[row]]++;

    for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
		for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
			h0[i][j] = t0[i][j] + t1[i][j];
}

void Table2SNPImpl2Phenotypes::fill(const UBYTE *encoded_genotype_pairs,
		const UBYTE *encoded_phenotypes, int size)
{

	int (*m)[NB_GENOTYPES_WITH_ZZ][NB_GENOTYPES_WITH_ZZ] = _matrix;
	int* T = &m[0][0][0];
	int (*t0)[NB_GENOTYPES_WITH_ZZ] = m[0];
	int (*t1)[NB_GENOTYPES_WITH_ZZ] = m[1];
	int (*h0)[NB_GENOTYPES_WITH_ZZ] = m[_H0_PHENOTYPE];


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

ostream & Table2SNPImpl2Phenotypes::printOn(std::ostream & o) const {
	for (int k = 0; k < _NB_PHENOTYPES_WITH_H0; ++k) {
		o << "phenotype=" << k << "\n";
		for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i) {
				for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
					o << _matrix[k][i][j] << "\t\t";
				o << "\n";
		}
	}

	return o;
}

