/*
 * divariate_table.hpp
 *
 * Latest benchmarks
 *
 * == more than 2 phenotypes (Table2SNPImplGeneric)
 * 1000 samples, 3 phenotypes, 2000000 times
 *
 * 	Table1SNp : 2.78s Table2Snps_std: 4.21 Table2Snps_opt: 2.23s
 *
 * N.B: Generic Table1SNp could be optimized, but in practice it is not needed
 * since Table1SNp are compute from Table2SNps
 *
 *  == more than 2 phenotypes (Table2SNPImpl2Phenotypes)
 *
 * 1000 samples, 2 phenotypes, 2000000 times
 * Table1SNp:2.07 	Table2Snps_std:2.76 	Table2Snps_opt:2.12
 *
 *
 *  Created on: Jun 22, 2010
 *      Author: kforner
 */

#ifndef TABLE2SNP_HPP_
#define TABLE2SNP_HPP_

#include "misc.hpp"
#include "matrix3d.hpp"

class Table2SNPImplInterface;
class Table2SNPImplGeneric;
class Table2SNPImpl2Phenotypes;

class Table2SNP {
public: // ========== LIFECYCLE ==========
	Table2SNP(int nb_phenotypes);
	~Table2SNP();

public: // =========== ACCESSORS ==========
	int nb_phenotypes() const;
	int const& operator()(int i, int j, int k) const;
	int& operator()(int i, int j, int k);

public: // =========== METHODS ==========

	void fill(const UBYTE* genotypes1, const UBYTE* genotypes2, const UBYTE* phenotypes, int size);
	void fill(vector<UBYTE> const& genotypes1, vector<UBYTE> const& genotypes2, vector<UBYTE> const& phenotypes);

	void fill(vector<UBYTE> const& encoded_genotype_pairs, vector<UBYTE> const& encoded_phenotypes);

	friend std::ostream& operator<< (std::ostream& o,  const Table2SNP& t);

protected:
	Table2SNPImplInterface& getTableImplementation();
	Table2SNPImplInterface const & getTableImplementation() const;

private: // ========== INSTANCE DATA ==========
	Table2SNPImplInterface* _table;
};


class Table2SNPImplInterface {
public: // ========== LIFECYCLE ==========
	virtual ~Table2SNPImplInterface() {};
public: // =========== ACCESSORS ==========
	virtual int nb_phenotypes() const = 0;
	virtual int const& operator()(int i, int j, int k) const = 0;
	virtual int& operator()(int i, int j, int k) = 0;
public: // =========== METHODS ==========
	virtual void fill(const UBYTE* genotypes1, const UBYTE* genotypes2, const UBYTE* phenotypes, int size) = 0;
	virtual void fill(const UBYTE* encoded_genotype_pairs, const UBYTE* encoded_phenotypes, int size) = 0;
	virtual ostream& printOn(std::ostream& o) const = 0;
};

class Table2SNPImplGeneric : public Table2SNPImplInterface{
public: // ========== LIFECYCLE ==========
	Table2SNPImplGeneric(int nb_phenotypes);
	virtual ~Table2SNPImplGeneric();
public: // =========== ACCESSORS ==========
	virtual int nb_phenotypes() const;
	virtual int const& operator()(int i, int j, int k) const { return _matrix(i,j,k); }
	virtual int& operator()(int i, int j, int k) { return _matrix(i,j,k); }
public: // =========== METHODS ==========
	void fill_for_H0();
	virtual void fill(const UBYTE* genotypes1, const UBYTE* genotypes2, const UBYTE* phenotypes, int size);
	virtual void fill(const UBYTE* encoded_genotype_pairs, const UBYTE* encoded_phenotypes, int size);
	virtual ostream& printOn(std::ostream& o) const;

private: // ========== INSTANCE DATA ==========
	int _nb_true_phenotypes;
	int _nb_phenotypes_with_h0;
	int _h0_phenotype;
	Matrix3D<int> _matrix;
};


class Table2SNPImpl2Phenotypes : public Table2SNPImplInterface {
public: // ========== LIFECYCLE ==========
	Table2SNPImpl2Phenotypes() {}
	// Based on the Law Of The Big Three:
	virtual ~Table2SNPImpl2Phenotypes() {}

public: // =========== ACCESSORS ==========
	virtual int nb_phenotypes() const;
	virtual int const& operator()(int i, int j, int k) const { return _matrix[i][j][k]; }
	virtual int& operator()(int i, int j, int k) { return _matrix[i][j][k]; }
//	int* getH0Row() { return _matrix[_H0_PHENOTYPE]; }
//	int const* getH0Row() const { return _matrix[_H0_PHENOTYPE]; }

public: // =========== METHODS ==========
	virtual void fill(const UBYTE* genotypes1, const UBYTE* genotypes2, const UBYTE* phenotypes, int size);
	virtual void fill(const UBYTE* encoded_genotype_pairs, const UBYTE* phenotypes, int size);
	virtual ostream& printOn(std::ostream& o) const;

private: // ========== CONSTANT ==========
	static const int _H0_PHENOTYPE = 2;
	static const int _NB_PHENOTYPES_WITH_H0 = 3;
	static const int _NB_PHENOTYPES = 2;
private: // ========== INSTANCE DATA ==========
	int _matrix[_NB_PHENOTYPES_WITH_H0][NB_GENOTYPES_WITH_ZZ][NB_GENOTYPES_WITH_ZZ];
};


// ========== Table2SNP inline IMPLEMENTATION ==========

inline int Table2SNP::nb_phenotypes() const {
	return getTableImplementation().nb_phenotypes();
}

inline Table2SNPImplInterface&
Table2SNP::getTableImplementation() { return *_table; }

inline void Table2SNP::fill(const UBYTE* genotypes1, const UBYTE* genotypes2, const UBYTE* phenotypes, int size)
{
	getTableImplementation().fill(genotypes1, genotypes2, phenotypes, size );
}

inline void Table2SNP::fill(vector<UBYTE> const& genotypes1, vector<UBYTE> const& genotypes2, vector<UBYTE> const& phenotypes)
{
	getTableImplementation().fill(&genotypes1[0], &genotypes2[0], &phenotypes[0], SIZE(phenotypes) );
}

inline void Table2SNP::fill(const vector<UBYTE> & encoded_genotype_pairs, const vector<UBYTE> & phenotypes)
{
	getTableImplementation().fill(&encoded_genotype_pairs[0], &phenotypes[0], SIZE(phenotypes) );
}



inline Table2SNPImplInterface const&
Table2SNP::getTableImplementation() const { return *_table; }

inline int& Table2SNP::operator ()(int i, int j, int k)
{
	return getTableImplementation()(i,j,k);
}

inline const int& Table2SNP::operator ()(int i, int j, int k) const
{
	return getTableImplementation()(i,j,k);
}

inline int Table2SNPImplGeneric::nb_phenotypes() const {
	return _nb_true_phenotypes;
}


inline int Table2SNPImpl2Phenotypes::nb_phenotypes() const {
	return _NB_PHENOTYPES;
}


#endif /* TABLE2SNP_HPP_ */
