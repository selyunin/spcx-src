/*
 * hybrid_automaton_visitor.cpp
 *
 *  Created on: Aug 31, 2009
 *      Author: frehse
 */

#include "core/hybrid_automata/hybrid_automaton_visitor.h"

#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/location.h"
#include "core/hybrid_automata/transition.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transform.h"
#include "utility/shared_ptr_output.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/hybrid_automata/automaton_id.h"

namespace hybrid_automata {

hybrid_automaton_visitor::~hybrid_automaton_visitor() {
}

void hybrid_automaton_visitor::prologue(hybrid_automaton& h) {
}

void hybrid_automaton_visitor::visit(hybrid_automaton& h, location& l, location_id l_id) {
}

void hybrid_automaton_visitor::visit(hybrid_automaton& h, transition& t, transition_id t_id) {
}

void hybrid_automaton_visitor::epilogue(hybrid_automaton& e) {
}

hybrid_automaton_visitor::visiting_order  hybrid_automaton_visitor::get_visiting_order() {
	return locations_before_transitions;
}

hybrid_automaton_visitor::network_visits hybrid_automaton_visitor::get_network_visits() {
	return all_components_and_subsets;
}


print_visitor::print_visitor(std::ostream& new_os) :
	os(new_os) {
}

print_visitor::~print_visitor() {
}

void print_visitor::prologue(hybrid_automaton& h) {
	os << "<hybrid_automaton>" << std::endl;
	os << "<name>" << h.get_name() << "</name>" << std::endl;
	os << "<declaration>" << std::endl;
	variable_id_set contrvars = h.get_variable_ids();
	variable_id_set inpvars = h.get_input_variables();
	set_difference_assign(contrvars, inpvars);
	for (variable_id_set::const_iterator it = contrvars.begin(); it != contrvars.end(); ++it) {
		os << "var " << variable(*it) << ";" << std::endl;
	}
	for (variable_id_set::const_iterator it = inpvars.begin(); it != inpvars.end(); ++it) {
		os << "input var " << variable(*it) << ";" << std::endl;
	}
	if (h.get_initial_states())
		os << "initial_states = " << h.get_initial_states() << ";" << std::endl;
	os << "</declaration>" << std::endl;
}

void print_visitor::epilogue(hybrid_automaton& h) {
	os << "</hybrid_automaton>" << std::endl;
}

void print_visitor::visit(hybrid_automaton& h, location& l, location_id l_id) {
	os << "<location id=\"id" << l_id << "\">" << std::endl;
	os << "<name>" << l.get_name() << "</name>" << std::endl;
	os << "<label kind =\"invariant\">";
	if(l.get_time_constraints().get_invariant())
	os << l.get_time_constraints().get_invariant();
	if(l.get_time_constraints().get_invariant() && l.get_time_constraints().get_dynamics())
	os << " & ";
	if(l.get_time_constraints().get_dynamics())
	os << l.get_time_constraints().get_dynamics();
	os << "</label>" << std::endl;
	os << "</location>" << std::endl;
}

void print_visitor::visit(hybrid_automaton& h, transition& t, transition_id t_id) {
	os << "<transition>" << std::endl;
	os << "<source ref=\"id" << t.get_source() << "\"/>" << std::endl;
	os << "<target ref=\"id" << t.get_target() << "\"/>" << std::endl;
	os << "<label kind =\"synchronisation\">" << named_label::get_name(t.get_label()) << "</label>"
			<< std::endl;
	os << "<label kind =\"assignment\">";
	if (t.get_jump_constraints().get_transform())
		os << t.get_jump_constraints().get_transform();
	os << "</label>" << std::endl;
	os << "<label kind =\"guard\">";
	if (t.get_jump_constraints().get_guard())
		os << t.get_jump_constraints().get_guard();
	os << "</label>" << std::endl;
	os << "</transition>" << std::endl;
}


print_visitor::network_visits print_visitor::get_network_visits() {
	return only_composition;
}

}

