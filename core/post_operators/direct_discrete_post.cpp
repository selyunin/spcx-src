/*
 * direct_discrete_post.cpp
 *
 *  Created on: Aug 31, 2009
 *      Author: frehse
 */

#include "core/post_operators/direct_discrete_post.h"

#include "core/continuous/continuous_set.h"
#include "core/continuous/continuous_set_operators.h"
#include "core/hybrid_automata/transition.h"
//#include "../abstract_framework/continuous/continuous_set_operators.h"

/** Forward declaration. */
namespace continuous {
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
continuous_set_ptr compute_intersection(const continuous_set_const_ptr& p1,
		const continuous_set_const_ptr& p2);
continuous_set_ptr compute_or_assign_intersection(continuous_set_ptr p1,
		const continuous_set_const_ptr& p2);
}


namespace hybrid_automata {

direct_discrete_post::~direct_discrete_post() {
}

direct_discrete_post* direct_discrete_post::clone() {
	return new direct_discrete_post(*this);
}

continuous::continuous_set_collection direct_discrete_post::post(const jump_constraints& trans,
		continuous::continuous_set::const_ptr source_inv,
		continuous::continuous_set::const_ptr target_inv,
		continuous::continuous_set::const_ptr cset) const {

	continuous::continuous_set_ptr ret_set;
	if (trans.get_guard()) // if there is a guard, intersect
	{
		// Assume constant_bound_dynamics
		ret_set = compute_intersection(cset, trans.get_guard());

		// Apply the transformation to intersection of cset with guard
		if (trans.get_transform())
			ret_set = compute_or_assign_transformation(ret_set,trans.get_transform());
	} else // no guard, go directly to transform
	{
		// Apply the transformation to cset
		if (trans.get_transform())
			ret_set = compute_transformation(cset,trans.get_transform());
		else {
			// no guard and no transform
			ret_set = continuous::continuous_set_ptr(cset->clone());
		}
	}

	IFLOGGER(DEBUG5) {
		std::stringstream ss;
		logger::copyfmt_to(ss);
		ss << ret_set;
		LOGGER(DEBUG5,"direct_discrete_post::post","result of map:\n"+ss.str());
	}

	// Intersect the result with the invariant of the target location
	if (target_inv) {
		IFLOGGER(DEBUG5) {
			std::stringstream ss;
			logger::copyfmt_to(ss);
			ss << target_inv;
			LOGGER(DEBUG5,"direct_discrete_post::post","adding target invariant "+ss.str());
		}

		ret_set = compute_or_assign_intersection(ret_set, target_inv);
	}

	if (ret_set->is_empty())
		LOGGER(DEBUG3,"direct_discrete_post::post","result of jump is empty");

	IFLOGGER(DEBUG5) {
		std::stringstream ss;
		logger::copyfmt_to(ss);
		ss << ret_set;
		LOGGER(DEBUG5,"direct_discrete_post::post","result of jump:\n"+ss.str());
	}

	return continuous::continuous_set_collection(ret_set);
}

}

