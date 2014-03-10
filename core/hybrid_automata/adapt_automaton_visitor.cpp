/*
 * adapt_automaton_visitor.cpp
 *
 *  Created on: Sep 11, 2009
 *      Author: frehse
 */

#include "core/hybrid_automata/adapt_automaton_visitor.h"

#include "utility/stl_helper_functions.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transforms.h"
#include "core/hybrid_automata/location.h"
#include "core/hybrid_automata/transition.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/symbolic_states/symbolic_state.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/predicates/dot_context_lookup.h"

namespace hybrid_automata {

adapt_discrete_set_visitor::adapt_discrete_set_visitor() {
	reset();
}

void adapt_discrete_set_visitor::reset() {
	success = true;
	my_discrete_set = discrete::discrete_set::ptr();
}

adapt_discrete_set_visitor::~adapt_discrete_set_visitor() {
}

bool adapt_discrete_set_visitor::get_success() {
	return success;
}

discrete::discrete_set::ptr adapt_discrete_set_visitor::get_discrete_set() {
	return my_discrete_set;
}

adapt_continuous_set_visitor::adapt_continuous_set_visitor() {
	reset();
}

void adapt_continuous_set_visitor::reset() {
	success = true;
	my_continuous_set = continuous::continuous_set::ptr();
}

adapt_continuous_set_visitor::~adapt_continuous_set_visitor() {
}

bool adapt_continuous_set_visitor::get_success() {
	return success;
}

continuous::continuous_set::ptr adapt_continuous_set_visitor::get_continuous_set() {
	return my_continuous_set;
}

adapt_dynamics_visitor::adapt_dynamics_visitor() :
	my_con_ad(adapt_continuous_set_visitor_ptr()) {
}

void adapt_dynamics_visitor::init(adapt_continuous_set_visitor_ptr con_ad) {
	my_con_ad = con_ad;
	my_invariant = continuous::continuous_set::ptr();
}

void adapt_dynamics_visitor::set_invariant(continuous::continuous_set::ptr inv) {
	my_invariant = inv;
}

void adapt_dynamics_visitor::reset() {
	success = true;
	my_dynamics = continuous::continuous_dynamics::ptr();
	if (my_con_ad)
		my_con_ad->reset();
}

adapt_dynamics_visitor::~adapt_dynamics_visitor() {
}

void adapt_dynamics_visitor::dispatch(
		const continuous::constant_bound_dynamics* d) {
	assert(my_con_ad);
	continuous::relation::const_ptr c = d->get_set();
	// convert the relation using continuous_set_adaptor
	my_con_ad->reset();
	c->accept(*my_con_ad);
	success = my_con_ad->get_success();
	continuous::relation::ptr new_c = my_con_ad->get_continuous_set();
	my_dynamics = continuous::constant_bound_dynamics::ptr(
			new continuous::constant_bound_dynamics(new_c));
}

void adapt_dynamics_visitor::dispatch(const continuous::relation_dynamics* d) {
	assert(my_con_ad);
	continuous::relation::const_ptr c = d->get_relation();
	// convert the relation using continuous_set_adaptor
	my_con_ad->reset();
	c->accept(*my_con_ad);
	success = my_con_ad->get_success();
	continuous::relation::ptr new_c = my_con_ad->get_continuous_set();
	my_dynamics = continuous::relation_dynamics::ptr(
			new continuous::relation_dynamics(new_c));
}

bool adapt_dynamics_visitor::get_success() {
	return success;
}

continuous::continuous_dynamics::ptr adapt_dynamics_visitor::get_dynamics() const {
	return my_dynamics;
}

continuous::continuous_set::ptr adapt_dynamics_visitor::get_invariant() const {
	return my_invariant;
}

adapt_transform_visitor::adapt_transform_visitor() :
	my_con_ad(adapt_continuous_set_visitor_ptr()) {
}

void adapt_transform_visitor::init(adapt_continuous_set_visitor_ptr con_ad) {
	my_con_ad = con_ad;
}

adapt_transform_visitor* adapt_transform_visitor::clone() {
	return new adapt_transform_visitor(*this);
}

void adapt_transform_visitor::reset() {
	success = true;
	my_transform = continuous::continuous_set_transform::ptr();
	if (my_con_ad)
		my_con_ad->reset();
}

adapt_transform_visitor::~adapt_transform_visitor() {
}

void adapt_transform_visitor::dispatch(const continuous::relation_transform* t) {
	assert(my_con_ad);
	continuous::relation::const_ptr c = t->get_relation_const();
	// convert the relation using continuous_set_adaptor
	my_con_ad->reset();
	c->accept(*my_con_ad);
	success = my_con_ad->get_success();
	continuous::relation::ptr new_c = my_con_ad->get_continuous_set();
	my_transform = continuous::relation_transform::ptr(
			new continuous::relation_transform(new_c));
}

void adapt_transform_visitor::dispatch(
		const continuous::constant_bound_time_elapse_transform* t) {
	assert(my_con_ad);
	continuous::relation::const_ptr c = t->get_set();
	// convert the relation using continuous_set_adaptor
	my_con_ad->reset();
	c->accept(*my_con_ad);
	success = my_con_ad->get_success();
	continuous::relation::ptr new_c = my_con_ad->get_continuous_set();
	my_transform = continuous::constant_bound_time_elapse_transform::ptr(
			new continuous::constant_bound_time_elapse_transform(new_c));
}
void adapt_transform_visitor::dispatch(
		const continuous::intersection_transform* t) {
	assert(my_con_ad);
	continuous::relation::const_ptr c = t->get_set();
	// convert the relation using continuous_set_adaptor
	my_con_ad->reset();
	c->accept(*my_con_ad);
	success = my_con_ad->get_success();
	continuous::relation::ptr new_c = my_con_ad->get_continuous_set();
	if (dynamic_cast<const continuous::pre_intersection_transform*> (t))
		my_transform = continuous::intersection_transform::ptr(
				new continuous::pre_intersection_transform(new_c));
	else
		my_transform = continuous::intersection_transform::ptr(
				new continuous::post_intersection_transform(new_c));
}
void adapt_transform_visitor::dispatch(const continuous::sequence_transform* t) {
	continuous::sequence_transform::ptr new_seq =
			continuous::sequence_transform::ptr(
					new continuous::sequence_transform());
	// Copy the transformed elements of the sequence into a new sequence
	for (continuous::sequence_transform::const_iterator it = t->begin(); it
			!= t->end(); ++it) {
		// create a new adaptor here instead of using the same one
		adapt_transform_visitor_ptr t_clone(clone());
		t_clone->reset();
		(*it)->accept(*t_clone);
		success = success && t_clone->get_success();
		new_seq->push_back(t_clone->get_transform());
	}
	my_transform = new_seq;
}

void adapt_transform_visitor::dispatch(
		const continuous::reset_affine_transform<global_types::rational_type>* t) {
	// convert the transform into a relation
	assert(my_con_ad);
	//@todo continuous::relation::ptr c = t->get_relation();
	continuous::relation::ptr c = continuous::relation::ptr();

	// convert the relation using continuous_set_adaptor
	my_con_ad->reset();
	c->accept(*my_con_ad);
	success = my_con_ad->get_success();
	continuous::relation::ptr new_c = my_con_ad->get_continuous_set();
	my_transform = continuous::relation_transform::ptr(
			new continuous::relation_transform(new_c));
}

void adapt_transform_visitor::dispatch(
		const continuous::reset_affine_transform<global_types::float_type>* t) {
	// convert the transform into a relation
	assert(my_con_ad);
	//@todo continuous::relation::ptr c = t->get_relation();
	continuous::relation::ptr c = continuous::relation::ptr();

	// convert the relation using continuous_set_adaptor
	my_con_ad->reset();
	c->accept(*my_con_ad);
	success = my_con_ad->get_success();
	continuous::relation::ptr new_c = my_con_ad->get_continuous_set();
	my_transform = continuous::relation_transform::ptr(
			new continuous::relation_transform(new_c));
}

void adapt_transform_visitor::dispatch(
		const continuous::reset_function_transform* t) {
	// convert the transform into a relation
	assert(my_con_ad);
	//@todo continuous::relation::ptr c = t->get_relation();
	continuous::relation::ptr c = continuous::relation::ptr();

	// convert the relation using continuous_set_adaptor
	my_con_ad->reset();
	c->accept(*my_con_ad);
	success = my_con_ad->get_success();
	continuous::relation::ptr new_c = my_con_ad->get_continuous_set();
	my_transform = continuous::relation_transform::ptr(
			new continuous::relation_transform(new_c));
}

bool adapt_transform_visitor::get_success() {
	return success;
}

continuous::continuous_set_transform::ptr adapt_transform_visitor::get_transform() {
	return my_transform;
}

adapt_automaton_visitor::adapt_automaton_visitor() :
	my_con_ad(adapt_continuous_set_visitor_ptr()), my_dyn_ad(
			adapt_dynamics_visitor_ptr()), my_trans_ad(
			adapt_transform_visitor_ptr()) {
}

void adapt_automaton_visitor::init() {
}

void adapt_automaton_visitor::define(adapt_continuous_set_visitor_ptr con_ad,
		adapt_dynamics_visitor_ptr dyn_ad,
		adapt_transform_visitor_ptr trans_ad, coefficient_type bool_t,
		coefficient_type number_t) {
	my_con_ad = con_ad;
	my_dyn_ad = dyn_ad;
	my_trans_ad = trans_ad;
	my_bool_type = bool_t;
	my_number_type = number_t;
	reset();
}

void adapt_automaton_visitor::reset() {
	my_loc_count = 0;
	my_trans_count = 0;
	my_success = true;
}

adapt_automaton_visitor::~adapt_automaton_visitor() {
}

adapt_automaton_visitor* adapt_automaton_visitor::clone() const {
	return new adapt_automaton_visitor(*this);
}

void adapt_automaton_visitor::visit(hybrid_automaton& h, location& l,
		location_id l_id) {
	// copy the time constraints in order to have exact copies of
	// properties other than dynamics and invariant
	time_constraints tcons = l.get_time_constraints();

	std::string problem_zone_msg = "";
	try {
		bool inv_ok = true;
		bool flow_ok = true;
		continuous::continuous_set::ptr adapted_inv =
				continuous::continuous_set::ptr();
		my_con_ad->reset();
		problem_zone_msg = "invariant of ";
		if (tcons.get_invariant()) {
			tcons.get_invariant()->accept(*my_con_ad);
			inv_ok = my_con_ad->get_success();
			if (inv_ok) {
				adapted_inv = my_con_ad->get_continuous_set();
				tcons.set_invariant(adapted_inv);
			}
		}
		problem_zone_msg = "flow of ";
		if (tcons.get_dynamics()) {
			my_dyn_ad->reset(); // resets also con_ad
			// pass the invariant to the dynamics adapter
			// so it can modify it if needed (quantify algebraic variables etc.)
			my_dyn_ad->set_invariant(adapted_inv);
			tcons.get_dynamics()->accept(*my_dyn_ad);
			flow_ok = my_dyn_ad->get_success();
			if (flow_ok) {
				tcons.set_dynamics(my_dyn_ad->get_dynamics());
				tcons.set_invariant(my_dyn_ad->get_invariant());
			}
		}
		if (flow_ok && inv_ok) {
			l.set_time_constraints(tcons);
			++my_loc_count;
			successful_visit();
		} else {
			failed_visit();
			throw basic_exception("Adaptation failure of undetermined cause.");
		}
	} catch (std::exception& e) {
		std::string loc = l.get_name();
		throw basic_exception("Failed to adapt " + problem_zone_msg
				+ "location " + loc + ".", e);
	}
}

void adapt_automaton_visitor::visit(hybrid_automaton& h, transition& t,
		transition_id t_id) {
	// copy the jump constraints in order to have exact copies of
	// properties other than guard and transform

	std::string problem_zone_msg = "";
	try {
		jump_constraints jcons = t.get_jump_constraints();
		my_trans_ad->reset();
		problem_zone_msg = "assignment of ";
		if (jcons.get_transform())
			jcons.get_transform()->accept(*my_trans_ad);
		my_con_ad->reset();
		problem_zone_msg = "guard of ";
		if (jcons.get_guard())
			jcons.get_guard()->accept(*my_con_ad);
		if (my_trans_ad->get_success() && my_con_ad->get_success()) {
			jcons.set_transform(my_trans_ad->get_transform());
			jcons.set_guard(my_con_ad->get_continuous_set());
			t.set_jump_constraints(jcons);
			++my_trans_count;
			successful_visit();
		} else {
			failed_visit();
			throw basic_exception("Adaptation failure of undetermined cause.");
		}
	} catch (std::exception& e) {
		std::string sloc = h.get_location(t.get_source())->get_name();
		std::string tloc = h.get_location(t.get_target())->get_name();
		std::string tlabel = dot_context::context_free_name(
				named_label::get_name(t.get_label()));
		throw basic_exception("Failed to adapt " + problem_zone_msg
				+ "transition from location " + sloc + " to " + tloc
				+ " with label " + tlabel + " (id " + to_string(t_id) + ").", e);
	}
}

bool adapt_automaton_visitor::get_success() const {
	return my_success;
}

void adapt_automaton_visitor::set_bool_type(coefficient_type t) {
	my_bool_type = t;
}
void adapt_automaton_visitor::set_number_type(coefficient_type t) {
	my_number_type = t;
}
;

void adapt_automaton_visitor::successful_visit() {
	my_success = my_success && true;
}
void adapt_automaton_visitor::failed_visit() {
	my_success = my_success && false;
}

adapt_continuous_set_visitor_ptr adapt_automaton_visitor::get_adapt_continuous_set_visitor() {
	return my_con_ad;
}
adapt_dynamics_visitor_ptr adapt_automaton_visitor::get_adapt_dynamics_visitor() {
	return my_dyn_ad;
}
adapt_transform_visitor_ptr adapt_automaton_visitor::get_adapt_transform_visitor() {
	return my_trans_ad;
}

adapt_automaton_visitor::coefficient_type adapt_automaton_visitor::get_bool_type() const {
	return my_bool_type;
}

adapt_automaton_visitor::coefficient_type adapt_automaton_visitor::get_number_type() const {
	return my_number_type;
}

void adapt_automaton_visitor::epilogue(hybrid_automaton& h) {
	// Adapt the initial set
	symbolic_state_collection::const_ptr ini_states = h.get_initial_states();
	if (ini_states) {
		try {
			//		symbolic_state_collection::ptr sstates = ini_states->create();
			//		for (symbolic_state_collection::const_iterator it = ini_states->begin(); it
			//				!= ini_states->end(); ++it) {
			//			// adapt the continuous set to get a new set
			//			my_con_ad->reset();
			//			(*it)->get_continuous_set()->accept(*my_con_ad);
			//			if (my_con_ad->get_success()) {
			//				symbolic_state::ptr s_new = symbolic_state::ptr(new symbolic_state(
			//						(*it)->get_discrete_set(), my_con_ad->get_continuous_set()));
			//				sstates->add(s_new);
			//			} else {
			//				failed_visit();
			//				return;
			//			}
			//		}
			symbolic_state_collection::ptr sstates(ini_states->clone());
			bool success = sstates->accept(*my_con_ad);
			if (success)
				successful_visit();
			else
				failed_visit();
			//	std::cout << "adapted ini states" << sstates;
			h.set_initial_states(sstates);
			//std::cout << "adapted ini states" << h.get_initial_states();
		} catch (std::exception& e) {
			throw basic_exception("Failed to adapt initial states.", e);
		}
	}
}

}

