/*
 * group_manager.hpp
 *
 * Simple utility class to manage groups. Currently manages a maximum of 256 groups.
 *
 *  Created on: Sep 22, 2011
 *
 */

#ifndef GROUP_MANAGER_HPP_
#define GROUP_MANAGER_HPP_

#include "misc.hpp"

/**
 * BEWARE: groups values should be < 256
 */
class GroupManager {
public: // ===== LIFECYCLE
	/**
	 * see <reset>
	 */
	GroupManager(const vector<int>& groups);
	GroupManager() {}
	~GroupManager() {}
public: // ==== INSTANCE METHODS =====
	int nb_groups() const { return SIZE(_indices_by_group); }
	const vector<int>& group_indices(int group) const { return _indices_by_group[group]; }

	/**
	 * reset/reinit the group manager for this array of group ids.
	 *
	 * @param a vector of group id/number per sample
	 *
	 */
	void reset(const vector<int>& groups);

private:
//	static const int MAX_BYTE =  1<<8;
	vector< vector<int> > _indices_by_group;
};


#endif /* GROUP_MANAGER_HPP_ */
