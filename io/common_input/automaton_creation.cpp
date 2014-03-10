#include "io/common_input/automaton_creation.h"

#include <iostream>
#include "core/continuous/polyhedra/ppl_polyhedron/convert_to_ppl.h"
#include "core/continuous/convert_predicate.h"
#include "core/continuous/predicate_continuous_set_constructors.h"
#include "core/symbolic_states/symbolic_state_collection_stl_list.h"
#include "core/continuous/continuous_dynamics/continuous_dynamics.h"
#include "core/continuous/predicate_continuous_set.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transforms.h"
#include "core/discrete/singleton_set.h"
#include "core/symbolic_states/symbolic_state.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/location_constraint_set.h"
#include "core/hybrid_automata/explicit_transition.h"

#include "io/common_input/predicate_parser.h"
#include "io/common_input/predicate_tree_operations.h"
#include "io/common_input/bool_node_creation.h"
#include "core/predicates/develop_matrix_expression.h"
#include "core/predicates/node_print_visitor.h"

namespace parser {
namespace automaton_parser {

using namespace hybrid_automata;

void add_variable(hybrid_automata::hybrid_automaton::ptr aut, std::string name,
		bool is_input, bool is_const, unsigned int dim) {

	variable_id id;
	/**
	 * If variable is a vector, add all element in automaton
	 */
	if (dim <= 1) {
		/** Instantiate scalars as vectors of dimension 1 */
		id = variable::get_or_add_variable_id(name, 1);
		/** But add only the vector variable to the automaton. */
		aut->add_variable(id, is_input, is_const);
	} else {
		/** Instantiate the vector */
		id = variable::get_or_add_variable_id(name, dim);
		for (int i = 1; i <= dim; ++i) {
			id = variable::get_or_add_variable_id(name + "(" + int2string(i)
					+ ")");
			aut->add_variable(id, is_input, is_const);
		}
	}
}

void add_symbol(hybrid_automata::hybrid_automaton_ptr aut, const symbol& symb) {
	if (symb.my_symbol_type != symbol::CONST_VALUE) {
		if (symb.my_symbol_type == symbol::LABEL) {
			// std::string lab = boost::any_cast<std::string>(symb.my_value(0, 0));
			std::string lab = symb.my_name;
			aut->add_label(lab);
		} else {
			bool is_constant = (symb.my_dynamics_type == symbol::CONSTANT);
			std::string name = symb.my_name;
			bool controlled = symb.is_controlled;
			unsigned int dim1 = symb.dim1;
			unsigned int dim2 = symb.dim2;
			unsigned int var_dim = 0;
			if (dim1 > 1 || dim2 > 1)
				var_dim = dim1 * dim2;
			add_variable(aut, name, !controlled, is_constant, var_dim);
		}
	}
}

tree::node::ptr extend_with_controlled_variable_relations(const tree::node::ptr& p,
		const variable_id_set& con_vars, const variable_id_set& uncon_vars, const variable_id_set& const_vars) {

	tree::node::ptr new_node = p;

	// Add const constraints for all variables in con_vars
	variable_id_set here_const_vars = con_vars;
	if (new_node) {
		variable_id_set mod_vars = valuation_functions::get_primed_variable_ids(new_node, 1);
		//std::cout << mod_vars;
		//operator<< <variable_id>(std::cerr,mod_vars);
		// don't make variables const that are modified by p
		set_difference_assign(here_const_vars, mod_vars);
	}
	// add const vars (parameters)
	set_union_assign(here_const_vars,const_vars);

	//		print_node(std::cerr,new_node.get());
	//std::cerr << " is modifying ";
	//print_variable_id_set(std::cerr,mod_vars);
	//std::cerr << " add const ";
	//print_variable_id_set(std::cerr,const_vars);
	tree::node::ptr const_rel = valuation_functions::create_const_relation(here_const_vars);
	new_node = valuation_functions::boolean_and(new_node,const_rel);

	// If there is still no assignment, but there are uncontrollable variables,
	// we need to create a proper node
	if (!new_node && !uncon_vars.empty()) {
		// There are uncontrolled variables but the transform is a null node.
		// Since the null node models all variables stay constant,
		// we must replace it with a "true" node, which is interpreted
		// as arbitrary reset.
		// @todo This shouldn't be a parser function.
		new_node = predicate_parser::create_bool_const_node("true");
	}

	return new_node;
}

void initialize_state(hybrid_automata::hybrid_automaton::ptr CurrentAutomaton,
		std::string str) {
	//cout<<"\n\n Initial State: "<<str<<"\n";
	//throw std::runtime_error("Missing implementation initialize_state");

	try {
		if (str != "") {
			// Construct a symbolic state

			// 1) Continuous part
			tree::node::ptr p =
					predicate_parser::create_bool_const_node("true");
			continuous::continuous_set_ptr cset =
					continuous::construct_predicate_continuous_set(p);

			// 2) Discrete part
			hybrid_automata::location_constraint_set lcs(
					CurrentAutomaton->get_id(),
					CurrentAutomaton->get_location_id(str));
			discrete::discrete_set::ptr dset = discrete::discrete_set::ptr(
					new discrete::singleton_set(lcs));

			// 3) Create the symbolic state
			symbolic_state::ptr sstate = symbolic_state::ptr(
					new symbolic_state(dset, cset));

			// 4) Create new symbolic_state_collection and add state
			symbolic_state_collection::ptr scol;
			if (CurrentAutomaton->get_initial_states()) {
				scol = symbolic_state_collection::ptr(
						CurrentAutomaton->get_initial_states()->create());
			} else {
				scol = symbolic_state_collection::ptr(
						new symbolic_state_collection_stl_list());
			}
			scol->add(sstate);
			CurrentAutomaton->set_initial_states(scol);
		}
	} catch (std::exception& e) {
		throw basic_exception("Failed to define initial states of automaton "
				+ CurrentAutomaton->get_name() + ".", e);
	}
}

hybrid_automata::location::ptr create_location(std::string loc_name,
		tree::node::ptr invariant, tree::node::ptr flow, const variable_id_set& const_vars) {

	hybrid_automata::location::ptr loc = hybrid_automata::location::ptr();
	try {
		// check that invariant don't contain constraints on primed variables
		variable_id_set invariant_mod_vars = valuation_functions::get_primed_variable_ids(invariant, 1);
		if (!invariant_mod_vars.empty()) {
			std::stringstream s;
			logger::copyfmt_to(s);
			print_variable_id_set(s,invariant_mod_vars);
			throw basic_exception("The derivatives of variable(s) "+s.str()+" are being assigned a value inside an invariant. They should be in the flow.");
		}

		// check that flow don't contain constraints without primed variables
		std::pair<tree::node_ptr, tree::node_ptr> flow_parts=divide_tree(flow);
		if (!valuation_functions::get_variable_ids(flow_parts.first).empty()) {
			std::stringstream s;
			logger::copyfmt_to(s);
			//s << flow_parts.first;
			print_as_predicate(s,flow_parts.first);
			throw basic_exception("The following constraints should be put inside an invariant instead of a flow: \n"+s.str());
		}

		invariant = valuation_functions::develop_assignment_matrix_tree(
				invariant);
		flow = valuation_functions::develop_assignment_matrix_tree(flow);

		// add const constraints for parameters
		valuation_functions::add_const_dynamics_constraints(flow,const_vars);

		continuous::continuous_set_ptr inv =
				continuous::construct_predicate_continuous_set(invariant);
		continuous::continuous_dynamics::ptr dyn =
				continuous::construct_predicate_continuous_set_dynamics(flow);

		hybrid_automata::time_constraints tcons(inv, dyn);
		loc = hybrid_automata::location::ptr(new hybrid_automata::location(
				loc_name, tcons));
	} catch (std::exception& e) {
		throw basic_exception("Failed to create location " + loc_name + ".", e);
	}

	return loc;
}

hybrid_automata::location::ptr create_location_from_string(
		std::string loc_name, std::string invariant,
		const variable_id_set& const_vars,
		const symbol_table& symbol_table, const parse_policy& ppol) {

	hybrid_automata::location::ptr lptr;
	std::pair<tree::node::ptr, tree::node::ptr> two_tree;
	try {
		tree::node::ptr p = predicate_parser::parse_predicate(invariant,
				symbol_table, ppol);
		two_tree = divide_tree(p);
	} catch (std::exception& e) {
		throw basic_exception("Failed to parse invariant or flow of location "
				+ loc_name + ".", e);
	}

	lptr = create_location(loc_name, two_tree.first, two_tree.second, const_vars);
	return lptr;
}

hybrid_automata::location::ptr create_location_from_string(
		std::string loc_name, std::string invariant, std::string flow,
		const variable_id_set& const_vars,
		const symbol_table& symbol_table, const parse_policy& ppol) {

	hybrid_automata::location::ptr lptr;
	tree::node::ptr inv;
	tree::node::ptr flo;
	try {
		inv = predicate_parser::parse_predicate(invariant, symbol_table, ppol);
	} catch (std::exception& e) {
		throw basic_exception("Failed to parse invariant of location "
				+ loc_name + ".", e);
	}
	try {
		flo = predicate_parser::parse_predicate(flow, symbol_table, ppol);
	} catch (std::exception& e) {
		throw basic_exception("Failed to parse flow of location " + loc_name
				+ ".", e);
	}

	lptr = create_location(loc_name, inv, flo, const_vars);
	return lptr;
}

hybrid_automata::transition::ptr create_transition(
		hybrid_automata::location_id sloc, std::string label,
		tree::node::ptr guard, tree::node::ptr assignment,
		hybrid_automata::location_id tloc, const variable_id_set& contr_vars,
		const variable_id_set& uncontr_vars, const variable_id_set& const_vars) {

	hybrid_automata::transition::ptr t;
	try {
		// check that guards don't contain constraints on primed variables
		variable_id_set guard_mod_vars = valuation_functions::get_primed_variable_ids(guard, 1);
		if (!guard_mod_vars.empty()) {
			std::stringstream s;
			logger::copyfmt_to(s);
			print_variable_id_set(s,guard_mod_vars);
			throw basic_exception("The variables "+s.str()+" are being assigned a value inside a guard. Put assignments in the assignment field.");
		}

		// check that assignments don't contain constraints without primed variables
		std::pair<tree::node_ptr, tree::node_ptr> assignment_parts=divide_tree(assignment);
		if (!valuation_functions::get_variable_ids(assignment_parts.first).empty()) {
			std::stringstream s;
			logger::copyfmt_to(s);
			print_as_predicate(s,assignment_parts.first);
			throw basic_exception("The following constraints should be put inside a guard instead of an assignment:\n"+s.str());
		}

//using ::operator<<;
//std::cout << "guard: " << guard << std::endl;
//std::cout << "assign: " << assignment << std::endl;

		assignment = valuation_functions::develop_assignment_matrix_tree(
				assignment);
//std::cout << "dev assign: " << assignment << std::endl;
		guard = valuation_functions::develop_assignment_matrix_tree(guard);
//std::cout << "extending ass with contr: ";
//print_variable_id_set(std::cout,contr_vars);
//std::cout << ", uncontr: ";
//print_variable_id_set(std::cout,uncontr_vars);
//std::cout << std::endl;
		assignment = extend_with_controlled_variable_relations(assignment,
				contr_vars, uncontr_vars, const_vars);
//std::cout << "ext assign: " << assignment << std::endl;

		continuous::continuous_set::ptr cguard =
				continuous::construct_predicate_continuous_set(guard);
		continuous::continuous_set_transform_ptr transf =
				continuous::construct_predicate_continuous_set_transform(
						assignment);

		hybrid_automata::jump_constraints jcons =
				hybrid_automata::jump_constraints(cguard, transf);
		t = hybrid_automata::transition::ptr(
				new hybrid_automata::explicit_transition(sloc, label, jcons,
						tloc));

	} catch (std::exception& e) {
		std::stringstream s;
		logger::copyfmt_to(s);
		throw basic_exception("Failed create transition with label " + label
				+ ".", e);
	}

	return t;
}

hybrid_automata::transition::ptr create_transition_from_string(
		hybrid_automata::location_id sloc, std::string label,
		std::string guard, std::string assignment,
		hybrid_automata::location_id tloc, const variable_id_set& contr_vars,
		const variable_id_set& uncontr_vars, const variable_id_set& const_vars, const symbol_table& symbol_table,
		const parse_policy& ppol) {

	tree::node::ptr guard_tree;
	try {
		guard_tree = predicate_parser::parse_predicate(guard, symbol_table, ppol);
	} catch (std::exception& e) {
		throw basic_exception("Failed to parse transition guard.", e);
	}

	tree::node::ptr assignment_tree;
	try {
		assignment_tree = predicate_parser::parse_predicate(assignment,
				symbol_table, ppol);
	} catch (std::exception& e) {
		throw basic_exception("Failed to parse transition assignment.", e);
	}
	return create_transition(sloc, label, guard_tree, assignment_tree, tloc,
			contr_vars, uncontr_vars, const_vars);
}

hybrid_automata::location::ptr dummy_location(std::string loc_name,
		std::string invariant) {
	hybrid_automata::location::ptr locp = create_location_from_string(loc_name,
			invariant);
	//	hybrid_automata::adapt_dynamics_visitor v;
	//	locp->get_time_constraints().get_dynamics()->accept(v);
	//	if (v.get_success()) {
	//		locp->get_time_constraints().set_dynamics(v.get_dynamics());
	//	} else
	//		throw std::runtime_error("Cannot adapt dynamics\n");
	return locp;
}
}

}
