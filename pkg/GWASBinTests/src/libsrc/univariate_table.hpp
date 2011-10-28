/*
 * univariatetable.hpp
 *
 *  Created on: Jun 23, 2010
 *      Author: kforner
 */

#ifndef UNIVARIATE_TABLE_HPP_
#define UNIVARIATE_TABLE_HPP_

#include "matrix.hpp"
#include "divariate_table.hpp"

class UnivariateTable  {
public: // ========== LIFECYCLE ==========
	UnivariateTable(int nb_phenotypes);

public: // =========== ACCESSORS ==========
	int nb_phenotypes() const;
	int nb_phenotypes_with_h0() const;
	int h0_phenotype() const;

	int const& operator()(int i, int j) const;
	int& operator()(int i, int j);

	int const* operator[](int i) const;
	int* operator[](int i);

	int const* getRow(int phenotype) const;

	int* getH0Row();
	int const* getH0Row() const;

	Matrix<int>& getMatrix();
	const Matrix<int>& getMatrix() const;

public: // =========== METHODS ==========
	int total(bool withMissingData) const;

	void fillForFirstSnp(const DivariateTable& t2);
	void fillForSecondSnp(const DivariateTable& t2);

	void fill_for_H0();
	void fill(vector<UBYTE> const& genotypes, vector<UBYTE> const& phenotypes);
	void fill(const UBYTE* genotypes, const UBYTE* phenotypes, int size);

	int computeSumBy(int phenotype) const;

	friend std::ostream& operator<< (std::ostream& o,  const UnivariateTable& t);

	bool operator ==(UnivariateTable const& t) const;
	bool operator !=(UnivariateTable const& t) const;

	void assign(int const (*m)[NB_GENOTYPES_WITH_ZZ]);

protected:
	void fill_nphenotypes(const UBYTE* genotypes, const UBYTE* phenotypes, int size);
	void fill_2phenotypes(const UBYTE* genotypes, const UBYTE* phenotypes, int size);

private: // ========== CONSTANTS ==========
	static const int NB_COLS = NB_GENOTYPES_WITH_ZZ;

private: // ========== INSTANCE DATA ==========
	int _nb_true_phenotypes;
	Matrix<int> _matrix;
};

inline void UnivariateTable::fill(const UBYTE *genotypes, const UBYTE *phenotypes, int size)
{
	if ( nb_phenotypes() == 2)
		fill_2phenotypes(genotypes, phenotypes, size);
	else
		fill_nphenotypes(genotypes, phenotypes, size);
}

inline void UnivariateTable::fill(const vector<UBYTE> & genotypes, const vector<UBYTE> & phenotypes)
{
	assert( genotypes.size() != 0 && phenotypes.size() == genotypes.size() );
	fill(&genotypes[0], &phenotypes[0], (int)genotypes.size() );
}


inline bool UnivariateTable::operator !=(UnivariateTable const& t) const {
	return ! (*this == t);
}

inline bool UnivariateTable::operator ==(UnivariateTable const& t) const {
	return _matrix == t._matrix;
}

inline int& UnivariateTable::operator ()(int i, int j)
{
	return _matrix[i][j];
}

inline const int& UnivariateTable::operator ()(int i, int j) const
{
	return _matrix[i][j];
}

inline int const* UnivariateTable::operator[](int i) const {
	return _matrix[i];
}

inline int* UnivariateTable::operator[](int i)  {
	return _matrix[i];
}

inline int UnivariateTable::total(bool withMissingData) const {
	const int *h0 = getH0Row();
	int sum =  h0[HOMOZYGOTE1] + h0[HETEROZYGOTE] + h0[HOMOZYGOTE2];
	if ( withMissingData )
		sum += h0[MISSING];
	return sum;
}

inline int UnivariateTable::nb_phenotypes() const {
	return _nb_true_phenotypes;
}

inline int UnivariateTable::h0_phenotype() const {
	return _nb_true_phenotypes;
}

inline int UnivariateTable::nb_phenotypes_with_h0() const {
	return _nb_true_phenotypes + 1;
}

inline int *UnivariateTable::getH0Row()
{
	return _matrix[h0_phenotype()];
}

inline const int *UnivariateTable::getRow(int phenotype) const
{
	assert( phenotype < nb_phenotypes_with_h0() );
	return _matrix[phenotype];
}

inline const int *UnivariateTable::getH0Row() const
{
	return _matrix[h0_phenotype()];
}

inline Matrix<int>& UnivariateTable::getMatrix() {
	return _matrix;
}
inline const Matrix<int>& UnivariateTable::getMatrix() const {
	return _matrix;
}


#endif /* UNIVARIATE_TABLE_HPP_ */
