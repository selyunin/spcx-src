/*
 * location_constraint_set_constructors.h
 *
 *  Created on: Dec 29, 2009
 *      Author: frehse
 */

#ifndef LOCATION_CONSTRAINT_SET_CONSTRUCTORS_H_
#define LOCATION_CONSTRAINT_SET_CONSTRUCTORS_H_

#include "core/hybrid_automata/location_constraint_set.h"
#include "core/predicates/valuation_function_tree_nodes.h"

/** Forward declarations */
namespace discrete {
class discrete_set;
typedef boost::shared_ptr<discrete_set> discrete_set_ptr;
}

namespace hybrid_automata {

/** Constructs a location_constraint_set from a conjunctive predicate.
 *
 * Assumes that the predicate is a conjunction.
 */
hybrid_automata::location_constraint_set construct_location_constraint_set(
		const tree::node::ptr& p);

/** Constructs a discrete set from a conjunctive predicate.
 *
 * Assumes that the predicate is a conjunction.
 * Creates a singleton_set as intermediate representation.
 */
discrete::discrete_set_ptr construct_discrete_set_from_conjunction(
		const tree::node::ptr& p);

}

#endif /* LOCATION_CONSTRAINT_SET_CONSTRUCTORS_H_ */
