/*
 * tables.hpp
 *
 *  Created on: Jun 18, 2010
 *      Author: kforner
 */

#ifndef TABLES1SNP_HPP_
#define TABLES1SNP_HPP_

#include "misc.hpp"
#include "matrix.hpp"
#include "table2snp.hpp"

class Table1SNPImplInterface;
class Table1SNPImpl;
class Table1SNPImpl2Phenotypes;

class Table1SNP {
public: // ========== LIFECYCLE ==========
	Table1SNP(int nb_phenotypes);
	~Table1SNP();

public: // =========== ACCESSORS ==========
	int nb_phenotypes() const;
	int const& operator()(int i, int j) const;
	int& operator()(int i, int j);

public: // =========== METHODS ==========
	void fill(const UBYTE* genotypes, const UBYTE* phenotypes, int size);
	void fill(vector<UBYTE> const& genotypes, vector<UBYTE> const& phenotypes);

	void fillForFirstSnp(const Table2SNP& t2);
	void fillForSecondSnp(const Table2SNP& t2);

	bool operator ==(Table1SNP const& t) const;
	bool operator !=(Table1SNP const& t) const;

	friend std::ostream& operator<< (std::ostream& o,  const Table1SNP& t);

protected:
	Table1SNPImplInterface& getTableImplementation();
	Table1SNPImplInterface const & getTableImplementation() const;

private: // ========== INSTANCE DATA ==========
	Table1SNPImplInterface* _table;
};



class Table1SNPImplInterface {
public: // ========== LIFECYCLE ==========
	//Table1SNPImplInterface(int nb_phenotypes);
	// Based on the Law Of The Big Three:
	virtual ~Table1SNPImplInterface() {};

public: // =========== ACCESSORS ==========
	virtual int nb_phenotypes() const = 0;
	virtual int const& operator()(int i, int j) const = 0;
	virtual int& operator()(int i, int j) = 0;
	virtual int* getH0Row() = 0;
	virtual int const* getH0Row() const = 0;

public: // =========== METHODS ==========

	virtual void fill(const UBYTE* genotypes, const UBYTE* phenotypes, int size) = 0;
	virtual ostream& printOn(std::ostream& o) const = 0;
};

class Table1SNPImplGeneric : public Table1SNPImplInterface{
public: // ========== LIFECYCLE ==========
	Table1SNPImplGeneric(int nb_phenotypes);
	// Based on the Law Of The Big Three:
	virtual ~Table1SNPImplGeneric();

public: // =========== ACCESSORS ==========
	virtual int nb_phenotypes() const;
	int const& operator()(int i, int j) const { return _matrix[i][j]; }
	int& operator()(int i, int j) { return _matrix[i][j]; }
	virtual int* getH0Row();
	virtual int const* getH0Row() const;
//	int const*const* getCArray() const;
//	int** getCArray();
public: // =========== METHODS ==========

	void fill_for_H0();
	virtual void fill(const UBYTE* genotypes, const UBYTE* phenotypes, int size);
	virtual ostream& printOn(std::ostream& o) const;

private: // ========== INSTANCE DATA ==========
	int _nb_true_phenotypes;
	int _nb_phenotypes_with_h0;
	int _h0_phenotype;
	Matrix<int> _matrix;
};


class Table1SNPImpl2Phenotypes : public Table1SNPImplInterface {
public: // ========== LIFECYCLE ==========
	Table1SNPImpl2Phenotypes() {}
	// Based on the Law Of The Big Three:
	virtual ~Table1SNPImpl2Phenotypes() {}

public: // =========== ACCESSORS ==========
	virtual int nb_phenotypes() const;
	int const& operator()(int i, int j) const { return _matrix[i][j]; }
	int& operator()(int i, int j)  { return _matrix[i][j]; }
	virtual int* getH0Row() { return _matrix[_H0_PHENOTYPE]; }
	virtual int const* getH0Row() const { return _matrix[_H0_PHENOTYPE]; }

public: // =========== METHODS ==========
	virtual void fill(const UBYTE* genotypes, const UBYTE* phenotypes, int size);
	virtual ostream& printOn(std::ostream& o) const;

private: // ========== CONSTANT ==========
	static const int _H0_PHENOTYPE = 2;
	static const int _NB_PHENOTYPES_WITH_H0 = 3;
	static const int _NB_PHENOTYPES = 2;
private: // ========== INSTANCE DATA ==========
	int _matrix[_NB_PHENOTYPES_WITH_H0][NB_GENOTYPES_WITH_ZZ];
};

// ========== Table1SNP inline IMPLEMENTATION ==========

inline Table1SNPImplInterface&
Table1SNP::getTableImplementation() { return *_table; }

inline Table1SNPImplInterface const&
Table1SNP::getTableImplementation() const { return *_table; }

inline bool Table1SNP::operator !=(Table1SNP const& t) const {
	return ! (*this == t);
}

inline int Table1SNP::nb_phenotypes() const {
	return getTableImplementation().nb_phenotypes();
}

inline int& Table1SNP::operator ()(int i, int j)
{
	return getTableImplementation()(i,j);
}

inline const int& Table1SNP::operator ()(int i, int j) const
{
	return getTableImplementation()(i,j);
}

inline void Table1SNP::fill(const UBYTE *genotypes, const UBYTE *phenotypes, int size)
{
	getTableImplementation().fill(genotypes, phenotypes, size);
}

inline void Table1SNP::fill(const vector<UBYTE> & genotypes, const vector<UBYTE> & phenotypes)
{
	assert( genotypes.size() != 0 && phenotypes.size() == genotypes.size() );
	fill(&genotypes[0], &phenotypes[0], (int)genotypes.size() );
}


// ========== Table1SNPImplGeneric inline IMPLEMENTATION ==========

inline int Table1SNPImplGeneric::nb_phenotypes() const {
	return _nb_true_phenotypes;
}

inline int *Table1SNPImplGeneric::getH0Row()
{
	return _matrix[_h0_phenotype];
}

inline const int *Table1SNPImplGeneric::getH0Row() const
{
	return _matrix[_h0_phenotype];
}

inline int Table1SNPImpl2Phenotypes::nb_phenotypes() const {
	return _NB_PHENOTYPES;
}


#endif /* TABLES1SNP_HPP_ */
