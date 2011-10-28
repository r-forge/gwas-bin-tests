/*
 * univariatetable.cpp
 *
 *  Created on: Jun 23, 2010
 *      Author: kforner
 */

#include "univariate_table.hpp"

UnivariateTable::UnivariateTable(int nb_phenotypes) :
		_nb_true_phenotypes(nb_phenotypes),
		_matrix(nb_phenotypes+1, NB_GENOTYPES_WITH_ZZ)
{}

void UnivariateTable::assign(int const (*m)[NB_GENOTYPES_WITH_ZZ]) {
	for (int k = 0; k < _nb_true_phenotypes; ++k)
		for (int i = 0; i < NB_COLS; ++i)
			_matrix[k][i] = m[k][i];
	fill_for_H0();
}

void UnivariateTable::fill_nphenotypes(const UBYTE *genotypes, const UBYTE *phenotypes, int size)
{
	_matrix.fill(0);
	Matrix<int>& m = _matrix;
	for (int i = 0; i < size; ++i) {
		m[phenotypes[i]][genotypes[i]]++;
		//m( phenotypes[i], genotypes[row] )++;
	}
	fill_for_H0();
}


void UnivariateTable::fill_2phenotypes(const UBYTE *genotypes, const UBYTE *phenotypes, int size)
{
	_matrix.fill(0);

	// seems to help the compiler optimize better
	int (*m)[4] = (int (*)[4])_matrix.c_contiguous_array();

	for (int row = 0; row < size; ++row)
		m[phenotypes[row]][genotypes[row]]++;

	int* h0 = getH0Row();

	for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
		h0[i] = m[0][i] + m[1][i];
}


void UnivariateTable::fillForFirstSnp(const DivariateTable& t2) {
	const int nb = nb_phenotypes();
	assert( t2.nb_phenotypes() == nb );

	const int* row = t2._matrix.c_contiguous_array();
	int* m = _matrix.c_contiguous_array();
	const int* mend = m + (nb+1)*NB_COLS;
	while (m != mend) {
		*m++ = row[0] + row[1] + row[2] + row[3];
			row += 4;
	}
}

void UnivariateTable::fillForSecondSnp(const DivariateTable& t2) {
	const int nb = nb_phenotypes();
	assert( t2.nb_phenotypes() == nb );
	UnivariateTable& t1 = *this;
	int sum = 0;
	for (int k = 0; k < nb; ++k) {
		for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i) {
			sum = 0;
			for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
				sum += t2(k,j,i);
			t1(k,i) = sum;
		}
	}

	fill_for_H0();
}

std::ostream& operator<< (std::ostream& o,  const UnivariateTable& t) {

	return o << t._matrix;
}

void UnivariateTable::fill_for_H0() {
	const int nb_phenotypes = _nb_true_phenotypes;
	int* h0_row = getH0Row();
	for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
		h0_row[i] = 0;
	int* row;
	for (int j = 0; j < nb_phenotypes; ++j) {
		row = _matrix[j];
		for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
			h0_row[i] += row[i];
	}
}


