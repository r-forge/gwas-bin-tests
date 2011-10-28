	/*
 * shuffle.cpp
 *
 *  Created on: Jun 14, 2010
 *      Author: kforner
 */

#include "group_manager.hpp"

using std::multimap;
using std::pair;

GroupManager::GroupManager(const vector<int>& groups) {
	reset(groups);
}

void GroupManager::reset(const vector<int> & groups)
{
	_indices_by_group.resize(0);

	multimap<int, int> hi;
	map<int, int> hk;

	const int nb = SIZE(groups);
	for (int i = 0; i < nb; ++i) {
		hi.insert(pair<int,int>(groups[i], i));
		hk[groups[i]]++;
	}

	int nb_groups = SIZE(hk);
	_indices_by_group.resize(nb_groups);
	int group_id = 0;
	map<int,int>::iterator it;
	map<int,int>::iterator ie = hk.end();
	pair<multimap<int,int>::iterator,multimap<int,int>::iterator> ret;
	multimap<int,int>::iterator mit;
	for (it = hk.begin() ; it != ie; it++ ) {
		int group = it->first;
		vector<int>& indices = _indices_by_group[group_id++];
		ret = hi.equal_range(group);
		for (mit=ret.first; mit!=ret.second; ++mit)
			indices.push_back(mit->second);
	}
/*
	const int nb = SIZE(groups);

	int hash_of_groups[MAX_BYTE];
	//  === compute the nb of groups ===
	for(int i = 0; i < MAX_BYTE; ++i)
		hash_of_groups[i] = 0;
	for (int i = 0; i < nb; ++i)
		hash_of_groups[ groups[i] ]++;
	int nb_groups = 0;
	for (int i = 0; i < MAX_BYTE; ++i)
		if ( hash_of_groups[i] )
			nb_groups++;

	//  === compute the nb of indices of each group: the group sizes ====

	_indices_by_group.resize(nb_groups);

	int current_group = 0;
	//  TRICK: get the nb in each group, then REPLACE it in the hash by the group number
	for (int i = 0; i < MAX_BYTE; ++i)
		if ( hash_of_groups[i] )
			hash_of_groups[i] = current_group++;

	for (int i = 0; i < nb; ++i) {
		int group_number = hash_of_groups[ groups[i] ];
		_indices_by_group[group_number].push_back(i);
	}
*/
}


