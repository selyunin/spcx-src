/*
 * hybrid_automaton_wrapper.cpp
 *
 *  Created on: Sep 10, 2009
 *      Author: frehse
 */

#include "core/hybrid_automata/hybrid_automaton_wrapper.h"

#include "utility/stl_helper_functions.h"
#include "utility/basic_exception.h"
#include "core/hybrid_automata/automaton_cache.h"
#include "core/hybrid_automata/location_constraint_set.h"

namespace hybrid_automata {

using namespace discrete;
using namespace continuous;

hybrid_automaton_wrapper::ptr hybrid_automaton_wrapper::get_ptr() {
	return boost::static_pointer_cast<hybrid_automaton_wrapper>(hybrid_automaton::get_ptr());
}

hybrid_automaton_wrapper::const_ptr hybrid_automaton_wrapper::get_const_ptr() const {
	return boost::static_pointer_cast<const hybrid_automaton_wrapper>(
			hybrid_automaton::get_const_ptr());
}

hybrid_automaton_wrapper::hybrid_automaton_wrapper() :
	my_impl(hybrid_automaton::ptr()) {
}

void hybrid_automaton_wrapper::set_impl(hybrid_automaton::ptr himpl) {
	assert(himpl);

//	if (my_impl) {
//		// swap identities back
//		hybrid_automaton_cache::swap_identity(get_ptr(), my_impl);
//	}

	// adopt new implementation
	my_impl = himpl;
	// set_name(my_impl->get_name() + "_WRAP" + int2string(my_impl->get_id()));
	// my_name=my_impl->my_name + "_IMPL" + int2string(my_id);

	// swap identities with implementation
//	hybrid_automaton_cache::swap_identity(get_ptr(), my_impl);

	assert(my_impl==himpl);
}

hybrid_automaton_wrapper::~hybrid_automaton_wrapper() {
}

hybrid_automaton* hybrid_automaton_wrapper::create() const {
	return my_impl->create();
}

void hybrid_automaton_wrapper::add_variable(const variable_id& vid, bool is_input, bool is_const) {
	my_impl->add_variable(vid, is_input, is_const);
}
void hybrid_automaton_wrapper::add_variables(const variable_id_set& vars,
		const variable_id_set& inp_vars,const variable_id_set& const_vars) {
	my_impl->add_variables(vars, inp_vars, const_vars);
}
const variable_id_set& hybrid_automaton_wrapper::get_variable_ids() const {
	return my_impl->get_variable_ids();
}
const variable_id_set& hybrid_automaton_wrapper::get_input_variables() const {
	return my_impl->get_input_variables();
}
const variable_id_set& hybrid_automaton_wrapper::get_const_variables() const {
	return my_impl->get_const_variables();
}

const symbolic_state_collection_ptr& hybrid_automaton_wrapper::get_initial_states() const {
	return my_impl->get_initial_states();
}

void hybrid_automaton_wrapper::set_initial_states(const symbolic_state_collection_ptr& sstate_set) {
	my_impl->set_initial_states(sstate_set);
}

std::pair<hybrid_automaton_wrapper::transition_const_iterator,
		hybrid_automaton_wrapper::transition_const_iterator> hybrid_automaton_wrapper::get_outgoing_transitions(
		location_id l, label_id a) const {
	return my_impl->get_outgoing_transitions(l, a);
}

transition_ptr hybrid_automaton_wrapper::get_transition(const transition_id& id) const {
	return my_impl->get_transition(id);
}

location_ptr hybrid_automaton_wrapper::get_location(const location_id& id) const {
	return my_impl->get_location(id);
}

location_id hybrid_automaton_wrapper::get_location_id(std::string loc_name) const {
	return my_impl->get_location_id(loc_name);
}

std::pair<hybrid_automaton_wrapper::location_const_iterator,
		hybrid_automaton_wrapper::location_const_iterator> hybrid_automaton_wrapper::get_locations() const {
	return my_impl->get_locations();
}

location_id_set hybrid_automaton_wrapper::get_locations(const location_constraint_set& lcons) const {
	// shortcut if it's a single constraint that defines the top implementation
	location_constraint_set cons=lcons;
	cons.map(get_id(),my_impl->get_id());
	return my_impl->get_locations(cons);
}

bool hybrid_automaton_wrapper::canonicalize_location_constraint(const automaton_id& aut_id, const location_constraint& con,
		location_constraint_set& lcons) const {
	automaton_id id = aut_id;
	if (aut_id == get_id())
		id = my_impl->get_id();
	return my_impl->canonicalize_location_constraint(id,con,lcons);
}

transition_id hybrid_automaton_wrapper::add_transition(const transition_ptr& trans, bool check_emptiness) {
	return my_impl->add_transition(trans, check_emptiness);
}

location_id hybrid_automaton_wrapper::add_location(const location_ptr& loc) {
	return my_impl->add_location(loc);
}

const label_id_set& hybrid_automaton_wrapper::get_labels() const {
	return my_impl->get_labels();
}

void hybrid_automaton_wrapper::add_label(const label_id& lab) {
//	if (lab == named_label::silent_id()) {
//		std::stringstream ss;
//		named_label::print_named_label_cache(ss);
//		throw basic_exception("silent label with id "+to_string(lab)+" should not be added to alphabet. Current global labels:\n"+ss.str());
//	}

	my_impl->add_label(lab);
}

void hybrid_automaton_wrapper::accept(hybrid_automaton_visitor& v) {
	std::string impl_name = my_impl->get_name();
	my_impl->set_name(get_name());

	my_impl->accept(v);

	my_impl->set_name(impl_name);
}

void hybrid_automaton_wrapper::print(std::ostream& os) const {
//	std::string impl_name=my_impl->get_name();
//	my_impl->set_name(get_name());
	my_impl->print(os);
//	my_impl->set_name(impl_name);
}

}
