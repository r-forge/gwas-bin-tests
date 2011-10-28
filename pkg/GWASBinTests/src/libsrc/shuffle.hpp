/*
 * shuffle.hpp
 *
 *  Created on: Jun 14, 2010
 *      Author: kforner
 */

#ifndef SHUFFLE_HPP_
#define SHUFFLE_HPP_

#include "misc.hpp"
#include "random.hpp"
#include "group_manager.hpp"


template <class T, class R>
void knuth_shuffle(T* array, int length, R& rand)
{
	int r, i = length;
	T t;
	while(i)
	{
		r = rand.random_int_below(i--);
		t = array[i];
		array[i] = array[r];
		array[r] = t;
	}
}

template<class T, class R>
void knuth_shuffle_by_indices(T* array, const int* indices, int indices_size, R& rand) {
	int r, i = indices_size;
	T t;
	int ind1, ind2;
	while (i) {
		r = rand.random_int_below(i--);
		ind1 = indices[i];
		ind2 = indices[r];
		t = array[ind1];
		array[ind1] = array[ind2];
		array[ind2] = t;
	}
}

template <class T, class R>
class Shuffler {
public:
	Shuffler(int seed1 = time(NULL), int seed2 = 0) :_rand(seed1, seed2) {}
	~Shuffler() {}
	inline void shuffle( vector<T>& v) {knuth_shuffle<T,R>(&v[0], SIZE(v), _rand);}
	inline void shuffle(T* array, int length)
	{ knuth_shuffle<T,R>(array, length, _rand); }
	inline void reseed(int seed1, int seed2 = 0) { _rand.reseed(seed1, seed2); }
protected:
	R _rand;
};

typedef Shuffler<UBYTE, RandomSFMT> StdShuffler;


template <class T, class R>
class GroupShuffler {
public:
	GroupShuffler(const vector<int>& groups, int seed1 = time(NULL), int seed2 = 0)
		: _gm(groups), _rand(seed1, seed2) {}
	~GroupShuffler() {}
	void shuffle( vector<T>& v) { shuffle(&v[0], SIZE(v)); }
	void shuffle(T* array, int length) {
		const int nb_groups = _gm.nb_groups();
		for(int i = 0; i < nb_groups; ++i) {
			const vector<int>& g = _gm.group_indices(i);
			knuth_shuffle_by_indices<T,R>(array, &g[0], SIZE(g), _rand);
		}
	}
	inline void reseed(int seed) { _rand.reseed(seed); }

protected:
	GroupManager _gm;
	R _rand;
};

typedef GroupShuffler<UBYTE, RandomSFMT> StdGroupShuffler;


#endif /* SHUFFLE_HPP_ */
