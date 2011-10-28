/*
 * tables.cpp
 *
 *  Created on: Jun 18, 2010
 *      Author: kforner
 */

#include "table1snp.hpp"



void Table1SNPImplGeneric::fill(const UBYTE *genotypes, const UBYTE *phenotypes, int size)
{
	_matrix.fill(0);
	Matrix<int>& m = _matrix;
	for (int i = 0; i < size; ++i) {
		m[phenotypes[i]][genotypes[i]]++;
		//m( phenotypes[i], genotypes[row] )++;
	}
	fill_for_H0();
}

void Table1SNPImplGeneric::fill_for_H0() {
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

ostream & Table1SNPImplGeneric::printOn(std::ostream & o) const
{
	return o << _matrix;
}

ostream & Table1SNPImpl2Phenotypes::printOn(std::ostream & o) const {
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
			o << _matrix[i][j] << "\t\t";
		o << "\n";
	}
	for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
		o << "--\t\t";
	o << "\n";
	for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
		o << _matrix[_H0_PHENOTYPE][j] << "\t\t";
	o << "\n";

	return o;
}

Table1SNPImplGeneric::Table1SNPImplGeneric(int nb_phenotypes) :
		_nb_true_phenotypes(nb_phenotypes),
		_nb_phenotypes_with_h0(nb_phenotypes+1),
		_h0_phenotype(nb_phenotypes),
		_matrix(nb_phenotypes+1, NB_GENOTYPES_WITH_ZZ)
{}

Table1SNPImplGeneric::~Table1SNPImplGeneric() {}

void Table1SNPImpl2Phenotypes::fill(const UBYTE *genotypes, const UBYTE *phenotypes, int size)
{
	int (*m)[4] = _matrix;
	for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
		m[0][i] = m[1][i] = m[2][i] = 0;

	for (int row = 0; row < size; ++row)
		m[phenotypes[row]][genotypes[row]]++;

	for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i)
		m[_H0_PHENOTYPE][i] = m[0][i] + m[1][i];
}

Table1SNP::Table1SNP(int nb_phenotypes)
{
	assert(nb_phenotypes >= 2);
	if ( nb_phenotypes == 2)
		_table = new Table1SNPImpl2Phenotypes();
	else
		_table = new Table1SNPImplGeneric(nb_phenotypes);
}



std::ostream& operator<< (std::ostream& o,  const Table1SNP& t) {
	return t._table->printOn(o);
}

Table1SNP::~Table1SNP()
{
	delete _table;
}

bool Table1SNP::operator ==(Table1SNP const& t) const {
	const int nb = nb_phenotypes();
	if ( t.nb_phenotypes() != nb )
		return false;
	for (int k = 0; k < nb; ++k) {
		for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i) {
			if ( (*this)(k,i) != t(k,i) )
				return false;
		}
	}

	return true;
}

void Table1SNP::fillForFirstSnp(const Table2SNP& t2) {
	const int nb = nb_phenotypes();
	assert( t2.nb_phenotypes() == nb );

	Table1SNP& t1 = *this;

	int sum = 0;
	for (int k = 0; k < nb; ++k) {
		for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i) {
			sum = 0;
			for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
				sum += t2(k,i,j);
			t1(k,i) = sum;
		}
	}
}

void Table1SNP::fillForSecondSnp(const Table2SNP& t2) {
	const int nb = nb_phenotypes();
	assert( t2.nb_phenotypes() == nb );

	Table1SNP& t1 = *this;

	int sum = 0;
	for (int k = 0; k < nb; ++k) {
		for (int i = 0; i < NB_GENOTYPES_WITH_ZZ; ++i) {
			sum = 0;
			for (int j = 0; j < NB_GENOTYPES_WITH_ZZ; ++j)
				sum += t2(k,j,i);
			t1(k,i) = sum;
		}
	}
}

//Table1SNPImpl::Table1SNPImpl(int nb_phenotypes)
//	: _nb_true_phenotypes(nb_phenotypes),
//	  _nb_phenotypes_with_h0(nb_phenotypes+1),
//	  _h0_phenotype(nb_phenotypes)
//{
//	_matrix = new Column[_nb_phenotypes_with_h0];
//}
//
//Table1SNPImpl::~Table1SNPImpl()
//{
//	delete[] _matrix;
//}
//
