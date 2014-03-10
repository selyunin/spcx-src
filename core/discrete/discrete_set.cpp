/*
 * discrete_set.cpp
 *
 *  Created on: Sep 1, 2009
 *      Author: frehse
 */

#include "core/discrete/discrete_set.h"

namespace discrete {

discrete_set::~discrete_set() {
}

bool discrete_set::is_empty() const {
	if (begin() == end())
		return true;
	else {
		for (const_iterator it = begin(); it != end(); ++it) {
			if (!it->is_empty()) return false;
		}
	}
	return true;
}

bool discrete_set::is_disjoint_from(const object_type& loc) const {
	for (const_iterator it = begin(); it != end(); ++it) {
		if (!it->is_disjoint_from(loc)) return false;
	}
	return true;
}

bool discrete_set::is_disjoint_from(const discrete_set::const_ptr& ds) const {
	for (const_iterator it = ds->begin(); it != ds->end(); ++it) {
		if (!is_disjoint_from(*it)) return false;
	}
	return true;
}

bool discrete_set::contains(const discrete_set::const_ptr& ds) const {
	// @todo make more efficient
	for (const_iterator it = ds->begin(); it != ds->end(); ++it) {
		if (!contains(*it)) return false;
	}
	return true;
}

}

