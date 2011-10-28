/*
 * genabel.cpp
 *
 *  Created on: Jun 16, 2010
 *      Author: kforner
 */

#include "genabel.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include "misc.hpp"
#include <cassert>


const int GenabelSnp::offsets[4] = { 6, 4, 2, 0 };
const char* GenabelSnp::allele_codes[]  = {"12", "AB" , "AT", "AG", "AC" , "A-" , "TA" ,
		"TG" , "TC", "T-", "GA" , "GT" , "GC" , "G-" , "CA", "CT", "CG" ,
		"C-", "-A", "-T", "-G" , "-C" , "21" , "BA"};

// ========== SETTERS ==========

void GenabelSnp::setChromosome(int chromosome) {
	this->chromosome = chromosome;
}

void GenabelSnp::setCoding(UBYTE coding) {
	this->coding = coding;
}

void GenabelSnp::setEncoded_genotypes(vector<UBYTE> const& encoded_genotypes) {
	this->encoded_genotypes = encoded_genotypes;
}

void GenabelSnp::setName(const string& name) {
	this->name = name;
}

void GenabelSnp::setPosition(int position) {
	this->position = position;
}

void GenabelSnp::setStrand(UBYTE strand) {
	this->strand = strand;
}


ostream & operator <<(ostream & s, const GenabelSnp & t)
{
	s << "GenabelSnp: " << t.getName() << " at "
			<< t.getChromosome() << ":" << (t.getStrand() ? '+' : '-')
			<< t.getPosition() << " for " << t.getNbGenotypes() << " genotypes\n";
	return s;
}

void GenabelSnp::printInfo(ostream& s) const {
	vector<UBYTE> genotypes;
	decodeGenotypes(genotypes);
	s << *this;

	string alleles;
	foreach(UBYTE g, genotypes) {
		genotypeToAlleles(g, alleles);
		s << alleles << ' ';
	}

	s << endl;
}

void GenabelSnp::check() const
{
	int max_size = SIZE(encoded_genotypes)*4;
	if ( nb_genotypes <= 0 || position < 0 ||
			nb_genotypes > max_size || nb_genotypes < (max_size - 3) )
	{
		cerr << "check failed for SNP " << name << endl;
	}

}

void GenabelSnp::genotype_to_allele(UBYTE genotype, UBYTE coding, string& dest)  {
	dest.resize(2);
	const char* code = allele_codes[coding - 1];
	switch (genotype) {
		case 0:
			dest[0] = code[0];
			dest[1] = code[0];
			break;
		case 1:
			dest[0] = code[0];
			dest[1] = code[1];
			break;
		case 2:
			dest[0] = code[1];
			dest[1] = code[1];
			break;

		default:
			dest[0] = 'Z';
			dest[1] = 'Z';
			break;
	}
}

void GenabelSnp::readEncodedGenotypesFromStream(istream & in, int nb_samples)
{

	const int nbytes = (int) ceil(nb_samples / 4.0);
	this->nb_genotypes = nb_samples;
	encoded_genotypes.clear();
	encoded_genotypes.reserve(nbytes);

	char buffer[3];
	buffer[2] = '\0';
	for (int i = 0; i < nbytes; ++i) {
		in >> buffer[0] >> buffer[1];
		encoded_genotypes.push_back( string_to_int(buffer, 16) );
	}
}

void GenabelSnp::readEncodedGenotypesFromArray(const UBYTE* encoded, int snp_index, int nb_samples)
{
	const int nbytes = (int) ceil(nb_samples / 4.0);
	int start = nbytes * snp_index;
	this->nb_genotypes = nb_samples;
	encoded_genotypes.resize(nbytes);
	std::copy(encoded + start, encoded + start + nbytes, encoded_genotypes.begin());


}

void GenabelSnp::decodeGenotypes(vector<UBYTE> & decoded) const
{
	decode_genotypes(getEncoded_genotypes(), getNbGenotypes(), decoded);
}


void GenabelSnp::decode_genotypes(const vector<UBYTE> & encoded, int nb_genotypes, vector<UBYTE> & decoded)
{
	const int nb_full_bytes = nb_genotypes/4;
	const int left_over = nb_genotypes % 4;
	decoded.clear();
	decoded.reserve(nb_genotypes);

	int b;
	int gt = 0;
	for (int i = 0; i < nb_full_bytes; ++i) {
		b = encoded[i];
		for (int j = 0; j < 4; ++j) {
			gt = ( (b>>offsets[j]) & 3);
			if ( gt > 0 )
				gt--;
			else
				gt = MISSING;
			decoded.push_back(gt);
		}
	}

	if ( left_over ) {
		b = encoded[nb_full_bytes];
		for (int j = 0; j < left_over; ++j) {
			gt = ( (b>>offsets[j]) & 3);
			if ( gt > 0 )
				gt--;
			else
				gt = MISSING;
			decoded.push_back(gt);
		}
	}
	assert((int)decoded.size() == nb_genotypes);
}


GenabelSnpData::GenabelSnpData(const string & genofile)
{
	load_snp_data(genofile);
	init();
}

GenabelSnpData::GenabelSnpData(const VS& sample_names,
		const VI& positions,
		const VU& codings,
		const VS& snp_names,
		const VI& chromosomes,
		const VU& strands,
		const VU& encoded_genotypes)
	: _sample_names(sample_names)
{

//	int nb_samples = SIZE(sample_names);
//	// sample names
//	_sample_names.reserve( nb_samples );
//	for (int i = 0; i < nb_samples; ++i) {
//		_sample_names.push_back( sample_names[i] );
//	}
	int nb_snps = SIZE(snp_names);
	_snps.reserve(nb_snps);
	const UBYTE* encoded = &encoded_genotypes[0];
	for (int i = 0; i < nb_snps; ++i) {
		GenabelSnp snp;
		snp.setName(snp_names[i]);
		snp.setChromosome(chromosomes[i]);
		snp.setCoding(codings[i]);
		snp.setPosition(positions[i]);
		snp.setStrand(strands[i]);
		snp.readEncodedGenotypesFromArray(encoded,i, SIZE(sample_names) );
		_snps.push_back( snp );
	}
	init();
}

void GenabelSnpData::check() const
{
	if ( nb_snps() <= 0 || nb_samples() <= 0 ) {
		cerr << "check failed for GenabelSnpData\n";
	}
	int nb = nb_snps();
	for (int i = 0; i < nb; ++i) {
		getSnp(i).check();
	}
}

void GenabelSnpData::init()
{
	// sort the snps
	sort(ALL(_snps));
//	const int nb = nb_snps();
//	for (int i = 0; i < nb; ++i) {
//		const GenabelSnp& snp = getSnp(i);
//		cerr << i << " " << snp.getName() << " " << snp.getChromosome() << " " << snp.getPosition() << endl;
//	}
}

void GenabelSnpData::load_snp_data(const string& genofile)
{

	const string VERSION = "version";
	const char delimiter = ' ';

	fstream in(genofile.c_str());
	if ( in.fail() ) {
		cerr << "failed opening file [" << genofile << "]\n";
		die("error opening snp_data file");
	}


	string header;
	getline(in, header);

	// fetch  version from header
	size_t pos = header.find( VERSION );
	if ( pos == string::npos)
		die("load_snp_data: bad header, could not find version info");
	istringstream hi( header.substr(pos + VERSION.length() + 1) );

	double  version;
	hi >> version;

	if (version <= 0) {
		cerr << "Bad version " << version << endl;
		die("load_snp_data: Bad version");
	}

	// N.B: lines of strings end with a space
	// === samples IDS ===
	string line;
	getline(in, line);
	int nids = count(line.begin(), line.end(), delimiter);
	split_line_by_space(line, _sample_names);
	int nb_samples = _sample_names.size();

	if( nb_samples != nids )
		die("problem with nb_samples");


	// === marker names ===
	getline(in, line);
	vector<string> snp_names;
	int nb_snps = count(line.begin(), line.end(), delimiter);
	snp_names.reserve(nb_snps);
	int n = split_line_by_space(line, snp_names);
	assert(n == nb_snps);

	// == chromosomes ==
	getline(in, line);
	vector<int> chromosomes;// need numeric type, UBYTE does not work
	chromosomes.reserve(nb_snps);
	n = split_line_by_space(line, chromosomes);
	assert( n == nb_snps );

	/* == positions == */
	getline(in, line);

	// I use doubles because there are positions coded in scientific notation e.g. 1.5e+07 that are not
	// well read by streams
	vector<double> positions;
	positions.reserve(nb_snps);
	n = split_line_by_space(line, positions);
	assert( n == nb_snps );
//
	// == SNP codings ==
	getline(in, line);
	n = line.size();
	vector<UBYTE> snp_codings;
	snp_codings.reserve(nb_snps);
	const char* begin = line.c_str();
	for (int i = 0; i < n; i += 3) {
		snp_codings.push_back( string_to_int(begin + i, 16) );
	}

	// == SNP strands ==
	getline(in, line);
	vector<short> strands; // need numeric type, UBYTE does not work
	n = split_line_by_space(line, strands);
	assert( n == nb_snps);

	// number of bytes to encode the genotypes (1 genotypes is stored on 2 bits

	for (int i = 0; i < nb_snps; ++i) {
		GenabelSnp snp;
		snp.setName(snp_names[i]);
		snp.setChromosome(chromosomes[i]);
		snp.setCoding(snp_codings[i]);
		snp.setPosition((int)positions[i]);
		snp.setStrand(strands[i]);
		snp.readEncodedGenotypesFromStream(in, nb_samples);
		_snps.push_back( snp );

	//	cerr << i << " " << snp.getName() << " " << snp.getChromosome() << " " << snp.getPosition() << endl;
	}

}

GenabelPhenoData::GenabelPhenoData(const string & phenofile)
{
	load_pheno_data(phenofile);
}

GenabelPhenoData::GenabelPhenoData(
		const vector<string>& ids
		,const vector<UBYTE>& sexs
		,const vector<UBYTE>& pops
		,const vector<int>& groups)
	: _ids(ids), _sexs(sexs), _pops(pops), _groups(groups)
{}


int GenabelPhenoData::max_groups() const
{
	if ( _groups.empty() )
		return 0;

	int nb_pop = *max_element(_groups.begin(), _groups.end()) + 1;
	return nb_pop;
}

void GenabelPhenoData::check() const
{
	int nb = SIZE(_ids);
	if ( nb != SIZE(_sexs) || nb != SIZE(_pops) ) {
		cerr << "check failed for GenabelPhenoData\n";
	}

}

void GenabelPhenoData::load_pheno_data(const string & genofile)
{
	fstream in(genofile.c_str());

	// ========== header ==========
	string header;
	getline(in, header);
	const char* HEADERS[] = { "id", "sex", "pop", "gws"};

	vector<string> headers;
	split_line_by_space(header, headers);


	// must contain  id sex pop and maybe gws IN THAT ORDER
	int nbc = headers.size();
	if ( nbc < 3 )
		die("Bad Number of columns");
	const char qq = '"';
	for (int i = 0; i < nbc; ++i) {
		removeSurroundingCharacter(headers[i], qq);
		if( HEADERS[i] != headers[i] ) {
			cerr << "Bad header column %" << (i+1) << ", should be " << HEADERS[i] << endl;
			die("Bad header");
		}
	}

	string line;
	while( in >> line ) {
		removeSurroundingCharacter(line, '"');
		_ids.push_back(line);

		in >> line;
		_sexs.push_back( string_to_int(line) );

		in >> line;
		_pops.push_back(string_to_int(line) );

		if ( nbc > 3) {
			in >> line;
			_groups.push_back( string_to_int(line) );
		}

	}
	if ( nbc <= 3)
		_groups.resize(_ids.size(), 0);

	assert( _ids.size() == _pops.size() );
	assert( _ids.size() == _sexs.size() );
	assert( _ids.size() == _groups.size() );

//	max_pop = *max_element(_pops.begin(), _pops.end());
//	max_sex = *max_element(_sexs.begin(), _sexs.end());
//	max_gws = *max_element(_groups.begin(), _groups.end());

}

// reorder the _ids[] to match the order of snp_data.getSampleNames()
void GenabelPhenoData::reorderBy(const GenabelSnpData& snp_data) {
	const int nb_samples = snp_data.nb_samples();

	if ( nb_samples != this->nb_samples() )
		die("Bad Number of samples");

	const vector<string>& samples = snp_data.getSampleNames();
	for (int i = 0; i < nb_samples; ++i) {

		int j = std::find(ALL(samples), this->_ids[i]) - samples.begin();
		if ( j == nb_samples ) { // not found
			cerr << "ERROR: id " << this->_ids[i] << " is not in Snp Data";
			die("Bad data");
		}

		// exchange the #i with #j
		if (i != j) {
			std::swap(_ids[i], _ids[j]);
			std::swap(_sexs[i], _sexs[j]);
			std::swap(_pops[i], _pops[j]);
			std::swap(_groups[i], _groups[j]);
		}

	}


}










