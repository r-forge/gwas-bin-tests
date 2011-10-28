/*
 * DivariateTable.hpp
 *
 *  Created on: Jun 23, 2010
 *      Author: kforner
 */

#ifndef DIVARIATETABLE_HPP_
#define DIVARIATETABLE_HPP_

#include "misc.hpp"
#include "matrix3d.hpp"


class UnivariateTable;

class DivariateTable {
	friend class UnivariateTable;

public: // ========== LIFECYCLE ==========
	DivariateTable(int nb_phenotypes);

public: // =========== ACCESSORS ==========
	int nb_phenotypes() const;
	int nb_phenotypes_with_h0() const;
	int h0_phenotype() const;

	int const*const* getMatrix(int phenotype) const;
	Matrix3D<int>& getMatrix();
	const Matrix3D<int>& getMatrix() const;

	int const& operator()(int i, int j, int k) const;
	int& operator()(int i, int j, int k);
	int const*const* operator[](int i) const;
	int** operator[](int i) ;
public: // =========== METHODS ==========
	int total(bool withMissingData) const;

	void fill_for_H0();
	void fill(const UBYTE* genotypes1, const UBYTE* genotypes2, const UBYTE* phenotypes, int size);
	void fill(vector<UBYTE> const& genotypes1, vector<UBYTE> const& genotypes2, vector<UBYTE> const& phenotypes);

	void fill(const UBYTE* encoded_genotype_pairs, const UBYTE* encoded_phenotypes, int size);
	void fill(vector<UBYTE> const& encoded_genotype_pairs, vector<UBYTE> const& encoded_phenotypes);

	bool operator ==(DivariateTable const& t) const;
	bool operator !=(DivariateTable const& t) const;

	friend std::ostream& operator<< (std::ostream& o,  const DivariateTable& t);

protected:
	void fill_nphenotypes(const UBYTE *encoded_genotype_pairs, const UBYTE *encoded_phenotypes, int size);

	void fill_2phenotypes(const UBYTE *encoded_genotype_pairs, const UBYTE *encoded_phenotypes, int size);

private: // ========== CONSTANTS ==========
	static const int NB_COLS = NB_GENOTYPES_WITH_ZZ;

private: // ========== INSTANCE DATA ==========
	int _nb_true_phenotypes;
	int _nb_phenotypes_with_h0;
	int _h0_phenotype;
	Matrix3D<int> _matrix;
};



// ========== DivariateTable IMPLEMENTATION ==========

inline int DivariateTable::total(bool withMissingData) const {
	const int*const* h0 = getMatrix(h0_phenotype());
	int nb = withMissingData ? NB_GENOTYPES_WITH_ZZ : NB_GENOTYPES_NO_ZZ;
	int sum = 0;
	for (int i = 0; i < nb; ++i)
		for (int j = 0; j < nb; ++j)
			sum += h0[i][j];

	return sum;
}

inline int DivariateTable::nb_phenotypes() const {
	return _nb_true_phenotypes;
}

inline int DivariateTable::h0_phenotype() const {
	return _nb_true_phenotypes;
}

inline int DivariateTable::nb_phenotypes_with_h0() const {
	return _nb_true_phenotypes + 1;
}

inline int const*const* DivariateTable::getMatrix(int phenotype) const {
	return _matrix[phenotype];
}

inline Matrix3D<int>& DivariateTable::getMatrix() {
	return _matrix;
}
inline const Matrix3D<int>& DivariateTable::getMatrix() const {
	return _matrix;
}


inline int const*const* DivariateTable::operator[](int i) const {
	return _matrix[i];
}

inline int** DivariateTable::operator[](int i)  {
	return _matrix[i];
}

inline bool DivariateTable::operator !=(DivariateTable const& t) const {
	return ! (*this == t);
}

inline int& DivariateTable::operator ()(int i, int j, int k)
{
	return _matrix(i,j,k);
}

inline const int& DivariateTable::operator ()(int i, int j, int k) const
{
	return _matrix(i,j,k);
}

inline void DivariateTable::fill(vector<UBYTE> const& genotypes1, vector<UBYTE> const& genotypes2,
		vector<UBYTE> const& phenotypes) {
	fill( &genotypes1[0], &genotypes2[0], &phenotypes[0], SIZE(genotypes1) );
}

inline void DivariateTable::fill(vector<UBYTE> const& encoded_genotype_pairs,
		vector<UBYTE> const& encoded_phenotypes) {
	fill( &encoded_genotype_pairs[0], &encoded_phenotypes[0], SIZE(encoded_genotype_pairs) );
}

inline void DivariateTable::fill(const UBYTE* encoded_genotype_pairs,
		const UBYTE* encoded_phenotypes, int size) {
	if ( nb_phenotypes() == 2 )
		fill_2phenotypes(encoded_genotype_pairs, encoded_phenotypes, size);
	else
		fill_nphenotypes(encoded_genotype_pairs, encoded_phenotypes, size);
}



#endif /* DIVARIATETABLE_HPP_ */
