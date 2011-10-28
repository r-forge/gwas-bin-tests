/*
 * bins.hpp
 *
 *  Created on: Jun 17, 2010
 *      Author: kforner
 */

#ifndef BINS_HPP_
#define BINS_HPP_

#include "misc.hpp"
#include "genabel.hpp"

class Bins {
public:
	Bins(const char* filename, const GenabelSnpData& snp_data, bool excludeX = false);
	Bins(const vector<string>& chr_names, vector<int>& vstart,  vector<int>& vends,
			const GenabelSnpData& snp_data, bool excludeX);

	friend ostream& operator <<(ostream& s, const Bins& t);


	void check() const;
	// ========== ACCESSORS ==========
	int nb_bins() const;
	int nb_non_empty_bins() const;
	int nb_snps_in_bins() const;

protected:
	// ========== INSTANCE METHODS ==========
	void load_bins(const char* filename);
	void compute_non_empty_bins(const GenabelSnpData& snp_data, bool filterX);

public: // ========== CLASS CONSTANTS ==========
//	static const int CHRX = 23;
//	static const int CHRY = 24;

public: // ========== CLASS METHODS ==========
//	static int chromosome_name_to_int(const string& name);
public:
	int nb_non_empty;
	int _nb_snps_in_bins;
	vector<int> chrs, starts, ends, first_snp, last_snp, non_empty_bins;
};


#endif /* BINS_HPP_ */
