/*
 * genabel.hpp
 *
 *  Created on: Jun 16, 2010
 *      Author: kforner
 */

#ifndef GENABEL_HPP_
#define GENABEL_HPP_

#include "misc.hpp"

class GenabelSnp {
public:
	GenabelSnp() : nb_genotypes(0), position(-1), coding(0), name(), chromosome(0), strand(0) {}

	// ========== INSTANCE METHODS ==========

	void decodeGenotypes(vector<UBYTE>& decoded) const;
	void readEncodedGenotypesFromStream(istream& in, int nb_samples);
	void readEncodedGenotypesFromArray(const UBYTE* encoded, int snp_index, int nb_samples);
	friend ostream& operator <<(ostream& s, const GenabelSnp& t);
	void printInfo(ostream& s) const;

	void genotypeToAlleles(UBYTE genotype, string& alleles) const;

	void check() const;

	// important, for sort: sort by genomic position
	bool operator < (const GenabelSnp& snp) const;

	// ========== GETTERS ==========

	int getNbGenotypes() const;
    int getChromosome() const;
    int getPosition() const;
    UBYTE getStrand() const;
    UBYTE getCoding() const;
    vector<UBYTE> getEncoded_genotypes() const;
	string const& getName() const;

    // ========== SETTERS ==========
    void setChromosome(int chromosome);
	void setCoding(UBYTE coding);
	void setEncoded_genotypes(vector<UBYTE> const& encoded_genotypes);
	void setName(const string& name) ;
	void setPosition(int position);
	void setStrand(UBYTE strand);

public: // ========== STATIC METHODS ==========
	static void decode_genotypes(const vector<UBYTE>& encoded, int nb_genotypes, vector<UBYTE>& decoded);
	static void genotype_to_allele(UBYTE genotype, UBYTE coding, string& dest);
private: // ========== CLASS CONSTANTS ==========
	static const int offsets[4];
	static const char* allele_codes[];

private: // ========== INSTANCE DATA ==========
	int nb_genotypes;
	int position;
	UBYTE coding;
	string name;
	int chromosome;
	UBYTE strand;
	vector<UBYTE> encoded_genotypes;
};



class GenabelSnpData {
	typedef vector<string> VS;
	typedef vector<int> VI;
	typedef vector<UBYTE> VU;
public:
	explicit GenabelSnpData(const string& genofile);
	GenabelSnpData(const VS& sample_names,
			const VI& positions,
			const VU& codings,
			const VS& snp_names,
			const VI& chromosomes,
			const VU& strands,
			const VU& encoded_genotypes);

	void check() const;
	void init();

public: // ========== GETTERS ==========
	const vector<string>& getSampleNames() const;
	const GenabelSnp& getSnp(int i) const;
    int nb_snps() const;
	int nb_samples() const;

protected:
	void load_snp_data(const string& genofile);

private:
	vector<string> _sample_names;
	vector< GenabelSnp > _snps;
};

class GenabelPhenoData {
public:
	GenabelPhenoData(const string& phenofile);
	GenabelPhenoData(const vector<string>& ids
			,const vector<UBYTE>& sexs
			,const vector<UBYTE>& pops
			,const vector<int>& groups = vector<int>(0));

public: // ========== GETTERS ==========
    vector<int> const& getGroups() const;
    vector<string> const& getIds() const;
    vector<UBYTE> const& getPops() const;
    vector<UBYTE> const& getSexs() const;

    int nb_samples() const;

public:	// ========== INSTANCE METHODS ==========
	void check() const;

    /*
     * return 0 if no group defined, otherwise the strict maximum of all group values
     * such as all group values are < max_groups()
     */
    int max_groups() const;

	void reorderBy(const GenabelSnpData& snp_data);

public: // ========== CLASS CONSTANTS ==========
	static const int FEMALE = 0;
	static const int MALE = 1;


protected:
	void load_pheno_data(const string& genofile);

private:
	vector<string> _ids;
	vector<UBYTE> _sexs;
	vector<UBYTE> _pops;
	vector<int> _groups;
};

// ========== GenabelSnp IMPLEMENTATION ==========

// ========== GETTERS ==========

inline void GenabelSnp::genotypeToAlleles(UBYTE genotype, string& alleles) const {
	genotype_to_allele(genotype, getCoding(), alleles);
}

inline int GenabelSnp::getNbGenotypes() const { return nb_genotypes; }

inline int GenabelSnp::getChromosome() const { return chromosome; }

inline int GenabelSnp::getPosition() const { return position; }

inline UBYTE GenabelSnp::getStrand() const { return strand; }

inline UBYTE GenabelSnp::getCoding() const { return coding;  }

inline vector<UBYTE> GenabelSnp::getEncoded_genotypes() const {return encoded_genotypes;}

inline string const& GenabelSnp::getName() const { return name; }

inline bool GenabelSnp::operator < (const GenabelSnp& snp) const {
	int c1 = getChromosome();
	int c2 = snp.getChromosome();
	return c1 < c2 ? true : (  c1 > c2 ? false : getPosition() < snp.getPosition() );
}


// ========== GenabelSnpData IMPLEMENTATION ==========
inline const vector<string>& GenabelSnpData::getSampleNames() const { return _sample_names; }

inline const GenabelSnp& GenabelSnpData::getSnp(int i) const {
	return _snps[i]; }

inline int GenabelSnpData::nb_snps() const { return _snps.size(); }

inline int GenabelSnpData::nb_samples() const { return _sample_names.size(); }

// ========== GenabelPhenoData IMPLEMENTATION ==========

inline vector<int> const& GenabelPhenoData::getGroups() const { return _groups; }
inline vector<string> const& GenabelPhenoData::getIds() const { return _ids; }
inline vector<UBYTE> const& GenabelPhenoData::getPops() const { return _pops; }
inline vector<UBYTE> const& GenabelPhenoData::getSexs() const { return _sexs; }
inline int GenabelPhenoData::nb_samples() const { return getIds().size(); }

#endif /* GENABEL_HPP_ */
