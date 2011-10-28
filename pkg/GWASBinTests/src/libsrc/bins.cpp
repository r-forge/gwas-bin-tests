/*
 * bins.cpp
 *
 *  Created on: Jun 17, 2010
 *      Author: kforner
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include <cerrno>

#include"bins.hpp"
#include "chromosomes.hpp"


Bins::Bins(const char* filename, const GenabelSnpData& snp_data, bool excludeX) {
	load_bins(filename);
	compute_non_empty_bins(snp_data, excludeX);
}

Bins::Bins(const vector<string>& chr_names, vector<int>& vstart,  vector<int>& vends, const GenabelSnpData& snp_data, bool excludeX)
 : chrs(), starts(vstart), ends(vends)
{
	chrs.reserve(chr_names.size());
	Chromosomes chr_converter;
	foreach(string s, chr_names)
		chrs.push_back( chr_converter.convert(s) );
	compute_non_empty_bins(snp_data, excludeX);
}

int Bins::nb_bins() const { return chrs.size(); }
int Bins::nb_non_empty_bins() const { return nb_non_empty; }
int Bins::nb_snps_in_bins() const { return _nb_snps_in_bins; }

ostream& operator <<(ostream& s, const Bins& b) {
	s << "Bins: non empty bins = " << b.nb_non_empty_bins() << '/' << b.nb_bins()
			<< ", " << b.nb_snps_in_bins() << " snps binned"<< endl;
	return s;
}

//int Bins::chromosome_name_to_int(const string& name) {
//	char c = toupper(name[0]);
//
//	switch (c)
//	{
//	case 'X':
//		return CHRX;
//	case 'Y':
//		return CHRY;
//	default:
//		break;
//	}
//
//	return string_to_int(name);
//}

void Bins::check() const
{
	int nb = nb_bins();
	if ( nb <= 0 || nb != SIZE(chrs) ||  nb != SIZE(starts) || nb != SIZE(ends)
			|| nb != SIZE(first_snp) || nb != SIZE(last_snp)  )
	{
		cerr << "check failed for Bins\n";
	}
}

void Bins::compute_non_empty_bins(const GenabelSnpData& snp_data, bool excludeX)
{

	const int nb_bins = this->nb_bins();
	const int nb_snps = snp_data.nb_snps();

	non_empty_bins.reserve(nb_bins);

	first_snp.resize(nb_bins, -1);
	last_snp.resize(nb_bins, -1);

	int snp_chr, bin_chr, snp_pos, bin_start, bin_end;
	int ibin = 0;
	int isnp = 0;
	int nb_snps_binned = 0;
	int nb_ne = 0;
	while( ibin < nb_bins && isnp < nb_snps ) {
		const GenabelSnp& snp = snp_data.getSnp(isnp);
		snp_pos = snp.getPosition();
		snp_chr = snp.getChromosome();
		bin_chr = chrs[ibin];
		bin_start = starts[ibin];
		bin_end = ends[ibin];

		// snp has no real position
		if ( snp_chr <= 0 || snp_pos <= 0) {
			isnp++;
			continue;
		}

		if ( excludeX && bin_chr == Chromosomes::X) {
			ibin++;
			continue;
		}

		// ## snp is in bin
		if ( snp_chr == bin_chr && snp_pos >= bin_start && snp_pos <= bin_end) {
			if (first_snp[ibin] < 0) {
				first_snp[ibin] = isnp;
				non_empty_bins.push_back(ibin);
				nb_ne++;
			}

			last_snp[ibin] = isnp;

			nb_snps_binned++;
			isnp++;
		} else if (snp_chr < bin_chr || (snp_chr == bin_chr && snp_pos < bin_start)) {
			/* snp is before bin ==> go to next snp */
			isnp++;
		} else if (snp_chr > bin_chr || (snp_chr == bin_chr && snp_pos > bin_end)) {
			/*  snp is after bin => go to next bin */
			ibin++;

		} else {
			die("should not happen");
		}
	}

	this->nb_non_empty = nb_ne;
	this->_nb_snps_in_bins = nb_snps_binned;
}

void Bins::load_bins(const char* filename) {
	fstream in(filename);
	if ( ! in ) {
		cerr << "unable to open file " << filename << endl;
		die("File Error");
	}

	Chromosomes chr_converter;

	string line;
	getline(in, line); // skip header

	while( in >> line ) {
		chrs.push_back( chr_converter.convert(line) );

		in >> line;
		starts.push_back( lexical_cast<int>(line) );

		in >> line;
		ends.push_back( lexical_cast<int>(line) );
	}
}
