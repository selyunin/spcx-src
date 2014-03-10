/*
 * predicate_continuous_set_constructors.cpp
 *
 *  Created on: Dec 29, 2009
 *      Author: frehse
 */

#include "predicate_continuous_set_constructors.h"
#include "core/continuous/continuous_dynamics/continuous_dynamics.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transforms.h"
#include "core/continuous/predicate_continuous_set.h"

namespace continuous {

continuous_set_ptr construct_predicate_continuous_set(const tree::node::ptr& p) {
	continuous_set_ptr cset = continuous_set_ptr(new predicate_continuous_set(p));
	return cset;
}

continuous_set_transform::ptr construct_predicate_continuous_set_transform(
		const tree::node::ptr& p) {
	continuous_set::ptr rel = construct_predicate_continuous_set(p);
	continuous_set_transform::ptr transf = continuous_set_transform::ptr(
			new relation_transform(rel));
	return transf;
}

continuous_dynamics::ptr construct_predicate_continuous_set_dynamics(const tree::node::ptr& p) {
	continuous_set_ptr deriv = construct_predicate_continuous_set(p);
	continuous_dynamics::ptr dyn = continuous_dynamics::ptr(new relation_dynamics(deriv));
	return dyn;
}

}

