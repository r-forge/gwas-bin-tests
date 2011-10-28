/*
 * misc.hpp
 *
 *  Created on: Jun 14, 2010
 *      Author: kforner
 */

#ifndef MISC_HPP_
#define MISC_HPP_



#include <ctime>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <cerrno>
#include <stdexcept>
#include <cassert>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/algorithm.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/casts.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/assign/std/vector.hpp>

#include "types.hpp"

// R stuff
// make sure Rprintf and REprintf are available
#ifndef USING_R
#define Rprintf(x) printf(x)
#define REprintf(x) fprintf(stderr,x)
#else
#include <R_ext/Print.h>
#endif

using namespace boost::assign;// bring 'operator+=()' into scope
using namespace boost::lambda;
using boost::lexical_cast;
using boost::bad_lexical_cast;

#define foreach BOOST_FOREACH
#define reverse_foreach BOOST_REVERSE_FOREACH
#define ALL(c) c.begin(), c.end()
#define SIZE(c) int((c).size())

using std::cout;
using std::cerr;
using std::cin;
using std::string;
using std::vector;
using std::istringstream;
using std::ostringstream;
using std::fstream;
using std::ostream;
using std::istream;
using std::endl;
using std::runtime_error;
using std::map;

template<typename K, typename V>
bool isInMap(const map<K,V>& map, const K& key) {
	return map.find(key) != map.end();
}

void removeSurroundingCharacter(string& s, char c);

int string_to_int(const char* s, int base=10) throw();

int string_to_int(const string& s, int base=10) throw();

// N.B: for lines with space delimitor with extra space at the end
template<typename T>
int split_line_by_space(const string& line, vector<T>& into)
{
	istringstream in(line);
	T s;
	char c;
	int cnt = 0;
	while ( in >> s) {
		into.push_back(s);
		cnt++;

		// TODO: check if this is still needed to consume the trailing delimiter
		in.get(c);
//		assert( c == ' ');
	}
	return cnt;
}

//template <class T>
//inline void swap(T& a, T& b)
//{
//	T t = a;
//	a = b;
//	b = t;
//}

void die(const char* msg);
double timespecDiff(const struct timespec *timeA_p, const struct timespec *timeB_p);

ostream& operator <<(ostream& s, const vector<UBYTE>& v);
template<typename T>
ostream& operator <<(ostream& o, const vector<T>& v)
{
	foreach(const T& b, v)
			o << b << ", ";
	//o.seekp(-1, std::ios_base::end);
	return o << endl;
}

class Timer {
public:
	Timer() { start(); }
	~Timer() {}
	void gettime(struct timespec* t) { clock_gettime(CLOCK_REALTIME, t); }
	void start() { gettime(&_start); }
	void stop() { gettime(&_end); }

	double elapsed(const struct timespec *timeA_p, const struct timespec *timeB_p) const
	{ return timespecDiff(timeA_p, timeB_p); }

	double elapsed() const {
		return elapsed(&_end, &_start);
	}

	friend ostream& operator <<(ostream& s, const Timer& t);
private:
	struct timespec _start;
	struct timespec _end;
};


#endif /* MISC_HPP_ */
