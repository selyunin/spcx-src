#ifndef AUTOMATON_CREATION_H_
#define AUTOMATON_CREATION_H_

#include <string>
#include "boost/shared_ptr.hpp"
#include "math/vdom/variable.h"
#include "core/hybrid_automata/location_id.h"
#include "io/common_input/symbol_table.h"
#include "parse_policy.h"

/** Forward declarations */
namespace tree {
class node;
typedef boost::shared_ptr<node> node_ptr;
}
namespace hybrid_automata {
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
typedef boost::shared_ptr<const hybrid_automaton> hybrid_automaton_const_ptr;
class transition;
typedef boost::shared_ptr<transition> transition_ptr;
typedef boost::shared_ptr<const transition> transition_const_ptr;
class location;
typedef boost::shared_ptr<location> location_ptr;
typedef boost::shared_ptr<const location> location_const_ptr;
}
namespace parser {
class symbol;
class symbol_table;
}


namespace parser {
namespace automaton_parser {

/**
 * Add the variable to the hybrid automaton.
 */
void add_variable(hybrid_automata::hybrid_automaton_ptr aut, std::string name, bool is_input=false, bool is_const = false, unsigned int dim = 0);

/**
 * Add the symbol (variable or label) to the hybrid automaton.
 *
 * If the symbol is a constant value, nothing is done.
 */
void add_symbol(hybrid_automata::hybrid_automaton_ptr aut, const symbol& symb);


/** Extend an assignment (jump relation) to controlled/uncontrolled variables
 *
 * Const constraints (z'==z) are added for all controlled variables (con_vars) that are not modified.
 * If there are uncontrolled variables but no controlled variables, a boolean node
 * "true" is added. This is necessary to model the fact that uncontrolled
 * variables can change their values arbitrarily (recall that a null node means
 * no change in any variable).
 *
 * \pre p must not contain OR or NOT nodes.
 *
 * @attention For uncontrolled variables to work correctly, a null node must be
 * replaced by an assignment if there are uncontrolled variables.
 * Otherwise, the uncontrolled variables would remain constant instead
 * of being quantified.
 */
tree::node_ptr extend_with_controlled_variable_relations(const tree::node_ptr& p,
		const variable_id_set& con_vars = variable_id_set(),
		const variable_id_set& uncon_vars = variable_id_set());

/**
 * Set the initial states of the automaton to be the location with name str.
 */
void initialize_state(
		hybrid_automata::hybrid_automaton_ptr CurrentAutomaton, std::string str);


/**
 * Create a location with two tree
 * @param loc_name the name of the location
 * @param invariant
 * @param flow
 * Convert trees to predicate_continuous_set and continuous_dynamics, then create the location.
 */
hybrid_automata::location_ptr create_location(
		std::string loc_name, tree::node_ptr invariant, tree::node_ptr flow, const variable_id_set& const_vars=variable_id_set());

/**
 * Create a location.
 * \param string loc_name the name of the location
 * \param string invariant the invariant, with invariant & flow predicates to parse
 */
hybrid_automata::location_ptr create_location_from_string(
		std::string loc_name, std::string invariant, const variable_id_set& const_vars=variable_id_set(), const symbol_table& symbol_table = symbol_table(), const parse_policy& ppol = parse_policy());

/**
 * Create a location.
 * \param string loc_name the name of the location
 * \param string invariant the invariant, with invariant to parse
 * \param string flow, with flow to parse
 */
hybrid_automata::location_ptr create_location_from_string(
		std::string loc_name, std::string invariant, std::string flow, const variable_id_set& const_vars=variable_id_set(), const symbol_table& symbol_table = symbol_table(), const parse_policy& ppol = parse_policy());

/**
 * Convert assignments to relation form
 * Convert trees to predicate_continuous_set and relation_transform
 * Combine them and create jump_constraints, then create the transition
 * @param sloc the source of transition
 * @param tloc the target of transition
 */
hybrid_automata::transition_ptr create_transition(
		hybrid_automata::location_id sloc, std::string label,
		tree::node_ptr guard, tree::node_ptr assignment,
		hybrid_automata::location_id tloc,
		const variable_id_set& contr_vars=variable_id_set(), const variable_id_set& uncontr_vars=variable_id_set(), const variable_id_set& const_vars=variable_id_set());

/**
 * Create a transition.
 * \param hybrid_automata::hybrid_automaton::location_id sloc the source of transition
 * \param hybrid_automata::hybrid_automaton::location_id tloc the target of transition
 * \param string guard the guard to parse
 * \param string assignment the assignment to parse
 */
hybrid_automata::transition_ptr create_transition_from_string(
		hybrid_automata::location_id sloc, std::string label,
		std::string guard, std::string assignment,
		hybrid_automata::location_id tloc, const variable_id_set& contr_vars =
				variable_id_set(), const variable_id_set& uncontr_vars =
				variable_id_set(), const variable_id_set& const_vars =
						variable_id_set(), const symbol_table& symbol_table =
				symbol_table(), const parse_policy& ppol = parse_policy());


/** Used by sfm routines to create locationes with ode_affine dynamics. */
hybrid_automata::location_ptr dummy_location(std::string loc_name, std::string invariant);

}
}

#endif
