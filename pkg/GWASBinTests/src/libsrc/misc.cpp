/*
 * misc.cpp
 *
 *  Created on: Jun 14, 2010
 *      Author: kforner
 */

#include "misc.hpp"

ostream& operator <<(ostream& o, const vector<UBYTE>& v)
{
	foreach(UBYTE b, v)
			o << int(b) << ", ";
	return o << endl;
}



int string_to_int(const char* s, int base) throw() {
	errno = 0;
	int i = strtol(s, NULL, base);
	if (errno != 0) {
		runtime_error error("not an integer : " + string(s));
		throw error;
	}
	return i;
}

int string_to_int(const string& s, int base) throw() {
	return string_to_int(s.c_str(), base);
}


void removeSurroundingCharacter(string& s, char c) {
	int begin = 0, end = s.size();
	if (s[0] == c)
		begin++;
	if ( s[end-1] == c )
		end --;
	s = s.substr(begin, end - begin);
}

//void swap(void** a, void** b)
//{
//	void* t = *a;
//	*a = *b;
//	*b = t;
//}
//void swap_UBYTE(UBYTE* a, UBYTE* b)
//{
//	UBYTE t = *a;
//	*a = *b;
//	*b = t;
//}

void die(const char* msg)
{
	fprintf(stderr, "AAAAARGHHHH ======> Died: %s\n", msg);
	assert(0);
	runtime_error error(msg);
	throw error;
}

double timespecDiff(const struct timespec *timeA_p, const struct timespec *timeB_p)
{
  return (((timeA_p->tv_sec * 1000000000) + timeA_p->tv_nsec) -
           ((timeB_p->tv_sec * 1000000000) + timeB_p->tv_nsec))/1000000000.0;
}

ostream & operator <<(ostream & s, const Timer & t)
{
	s << "Elapsed time: " << std::setprecision(2) << t.elapsed() << " s\n";
	return s;
}


