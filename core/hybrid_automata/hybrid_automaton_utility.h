/*
 * hybrid_automaton_utility.h
 *
 *  Created on: Sep 21, 2009
 *      Author: frehse
 */

#ifndef HYBRID_AUTOMATON_UTILITY_H_
#define HYBRID_AUTOMATON_UTILITY_H_

#include "boost/shared_ptr.hpp"
#include "core/hybrid_automata/location_constraint_set.h"

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
typedef boost::shared_ptr<const hybrid_automaton> hybrid_automaton_const_ptr;
class symbolic_state_collection;
typedef boost::shared_ptr<symbolic_state_collection> symbolic_state_collection_ptr;
typedef boost::shared_ptr<const symbolic_state_collection> symbolic_state_collection_const_ptr;
class adapt_automaton_visitor;
}

namespace hybrid_automata {

/** Canonicalize the location constraints in lcons.
 *
 * See hybrid_automaton.h for more information on canonic constraints. */
location_constraint_set canonicalize_location_constraints(
		const location_constraint_set& lcons);

/** Canonicalize the location constraints in a symbolic_state_collection.
 *
 * See hybrid_automaton.h for more information on canonic constraints. */
bool canonicalize_location_constraints(symbolic_state_collection& s);

/** Adapt the continuous sets in sstates according to v. */
bool adapt(adapt_automaton_visitor& v, symbolic_state_collection& sstates);

/** Intersect a symbolic_state_collection with the invariants of the
 * automata it refers to.
 *
 * The symbolic states are interpreted via the automaton aut,
 * so normally aut is a network and the symbolic states refer to
 * the automata in the network.
 *
 * The time equivalent discrete sets in s are enumerated, and for each
 * one the continuous set is intersected with the corresponding invariant.
 * Since the same enumeration is done for time elapse of the initial states,
 * it seems that this is the best one can do.
 *
 * @note The easy way to implement this would be to obtain the invariants as a
 * symbolic_state_collection, and then simply intersect the two collections.
 * However, this might mean enumerating an extremely large set of discrete states,
 * since the discrete state space is exponential in the number of automata
 * in the network.
 */
symbolic_state_collection_ptr intersect_with_invariants(
		const symbolic_state_collection& s,
		const hybrid_automaton& aut);

/** Convert a location_id to canonic location constraints */
location_constraint_set canonicalize_location_constraint(
		const hybrid_automaton_const_ptr& aut, const location_id& loc_id);

}

#endif /* HYBRID_AUTOMATON_UTILITY_H_ */
