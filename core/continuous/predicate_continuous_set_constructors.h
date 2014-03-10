/*
 * predicate_continuous_set_constructors.h
 *
 *  Created on: Dec 29, 2009
 *      Author: frehse
 */

#ifndef PREDICATE_CONTINUOUS_SET_CONSTRUCTORS_H_
#define PREDICATE_CONTINUOUS_SET_CONSTRUCTORS_H_

#include "core/predicates/valuation_function_tree_nodes.h"

/** Forward declarations */
namespace continuous {
class continuous_set;
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
class continuous_set_transform;
typedef boost::shared_ptr<continuous_set_transform> continuous_set_transform_ptr;
class continuous_dynamics;
typedef boost::shared_ptr<continuous_dynamics> continuous_dynamics_ptr;
}

namespace continuous {

/** Converts a predicate to a continuous set.
 *
 * Creates a predicate_continuous_set as intermediate representation.
 */
continuous_set_ptr construct_predicate_continuous_set(const tree::node::ptr& p);

/** Converts a predicate to a continuous transform.
 *
 * Creates a relation_transform as intermediate representation, where the relation
 * is represented as a continuous_set.
 */
continuous_set_transform_ptr
construct_predicate_continuous_set_transform(const tree::node::ptr& p);

/** Converts a predicate to a continuous dynamics.
 *
 * Creates a relation_dynamics as intermediate representation, where the relation
 * is represented as a continuous_set.
 */
continuous::continuous_dynamics_ptr
		construct_predicate_continuous_set_dynamics(const tree::node::ptr& p);

}

#endif /* PREDICATE_CONTINUOUS_SET_CONSTRUCTORS_H_ */
