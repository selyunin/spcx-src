/*
 * composition_operators.h
 *
 *  Created on: Aug 26, 2009
 *      Author: frehse
 */

#ifndef COMPOSITION_OPERATORS_H_
#define COMPOSITION_OPERATORS_H_

#include "core/hybrid_automata/transition.h"
/** Forward declarations */
namespace hybrid_automata {
class location;
typedef boost::shared_ptr<location> location_ptr;
typedef boost::shared_ptr<const location> location_const_ptr;

class symbolic_state;
typedef boost::shared_ptr<symbolic_state> symbolic_state_ptr;
typedef boost::shared_ptr<const symbolic_state> symbolic_state_const_ptr;
class symbolic_state_collection;
typedef boost::shared_ptr<symbolic_state_collection> symbolic_state_collection_ptr;
typedef boost::shared_ptr<const symbolic_state_collection> symbolic_state_collection_const_ptr;

class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
typedef boost::shared_ptr<const hybrid_automaton> hybrid_automaton_const_ptr;
}


namespace hybrid_automata {

/** Compute the set of input variables of the composition of two automata h1 and h2. */
variable_id_set compute_input_variables(const hybrid_automaton_const_ptr& h1, const hybrid_automaton_const_ptr& h2);

/** Compute the set of const-dynamics variables of the composition of two automata h1 and h2. */
variable_id_set compute_const_variables(const hybrid_automaton_const_ptr& h1, const hybrid_automaton_const_ptr& h2);

/** Create a new location that results from the parallel composition of
 * locations *p1 and *p2 */
location_ptr compose_location(const location_const_ptr& p1, const location_const_ptr& p2);

/** Compute the new jump constraints for two synchronized transitions with
 * constraints c1 and c2. */
jump_constraints compose_jump_constraints(const jump_constraints& c1, const jump_constraints& c2);

/** Compute the new jump constraints for an autonomous (non-synchronized) transition, in which
 * the second automaton h2 does not participate.
 *
 * c1 is modified to ensure that the controlled (non-input) variables of h2 remain constant.
 * For non-I/O semantics, don't force constness on variables that are controlled by h1 as well.
 * In an I/O framework, variables are controlled by exactly one automaton, so
 * this rule doesn't interfere with I/O. */
jump_constraints compose_autonomous_jump_constraints(
		const jump_constraints& c1, const hybrid_automaton_const_ptr& h1,
		const hybrid_automaton_const_ptr& h2);

/** Compute the composition of the symbolic states, i.e., the symbolic states
 * as interpreted on the composed automaton.
 *
 * If s1 is null, then s2 is returned and vice versa.
 *  */
symbolic_state_collection_ptr compose_symbolic_states(const symbolic_state_collection_const_ptr& s1,
		const symbolic_state_collection_const_ptr& s2);

}

#endif /* COMPOSITION_OPERATORS_H_ */
