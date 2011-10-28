/*
 * chromosomes.hpp: Mainly converting chromosome names to integer codes
 *
 *  Created on: Sep 3, 2010
 *      Author: kforner
 */

#ifndef CHROMOSOMES_HPP_
#define CHROMOSOMES_HPP_

#include "misc.hpp"

class Chromosomes {
public:
	Chromosomes();

public: // ========== INSTANCE METHODS =====
	int convert(const string& s);

public: // ========== CLASS METHODS =====

	bool isAutosomal(int chr) { return chr < X; }
	// return ERROR(0) on error

	int convert(int chr) {
		if ( chr > 0 && chr <= MT )
			return chr;
		return CHR_ERROR;
	}

public: // ====== CONSTANTS ======
	static const int X = 23;
	static const int Y = 24;
	static const int XY = 25;
	static const int MT = 26;
	static const int CHR_ERROR = 0;

private: // ===== INSTANCE DATA =====
	map<string, int> _names_to_chr;
};




#endif /* CHROMOSOMES_HPP_ */
