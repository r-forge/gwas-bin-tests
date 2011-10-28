
/*
 * chromosomes.cpp
 *
 *  Created on: Sep 3, 2010
 *      Author: kforner
 */

#include "chromosomes.hpp"

Chromosomes::Chromosomes()
{
	_names_to_chr["X"] = X;
	_names_to_chr["x"] = X;

	_names_to_chr["Y"] = Y;
	_names_to_chr["y"] = Y;

	_names_to_chr["XY"] = XY;
	_names_to_chr["xy"] = XY;


	_names_to_chr["MT"] = MT;
	_names_to_chr["mt"] = MT;

	for(int i = 0; i < 23; i++) {
		ostringstream out;
		out << i;
		_names_to_chr[out.str()] = i;

	}

}


int Chromosomes::convert(const string & s)
{
	if ( ! isInMap(_names_to_chr, s))
		return CHR_ERROR;
	return _names_to_chr[s];
//	int l = s.length();
//	if ( l <= 0 || l > 2)
//		return ERROR;
//	char s1 = toupper(s[0]);
//
//	if ( isdigit(s1) )
//		return string_to_int(s);
//
//	char s2 = toupper(s[1]);
}


