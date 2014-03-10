/*
 * ppl_adaptors.cpp
 *
 *  Created on: Sep 11, 2009
 *      Author: frehse
 */

#include "core/hybrid_automata/adaptors/ppl_adaptors.h"

#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"
#include "core/continuous/polyhedra/ppl_polyhedron/continuous_set_PPL_NNC.h"
#include "core/continuous/polyhedra/ppl_polyhedron/convert_to_ppl.h"

#include "core/continuous/continuous_dynamics/continuous_dynamics.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transforms.h"
#include "core/continuous/predicate_continuous_set.h"

#include "math/ode_solving/traj_simu/continuous_set_simulation.h"

#include "core/hybrid_automata/hybrid_automaton.h"

namespace hybrid_automata {

void continuous_set_to_PPL_adaptor::dispatch(const continuous::support_function_provider* c) {
	my_continuous_set = continuous::continuous_set::ptr();
	success = false;
}

void continuous_set_to_PPL_adaptor::dispatch(const continuous::constr_polyhedron<global_types::rational_type>* c) {
	ppl_polyhedron::continuous_set_PPL_NNC::ptr new_c =
			ppl_polyhedron::continuous_set_PPL_NNC::ptr(
					new ppl_polyhedron::continuous_set_PPL_NNC());
	new_c->add_constraints(c->get_constraints());
	my_continuous_set=new_c;
	success = true;
}

void continuous_set_to_PPL_adaptor::dispatch(const continuous::constr_polyhedron<global_types::float_type>* c) {
	my_continuous_set = continuous::continuous_set::ptr();
	success = false;
}

void continuous_set_to_PPL_adaptor::dispatch(const ppl_polyhedron::continuous_set_PPL_NNC* c) {
	// Nothing to be done
	// @todo this const-cast is dirty, but it saves us from copying the object
	ppl_polyhedron::continuous_set_PPL_NNC* nc=const_cast<ppl_polyhedron::continuous_set_PPL_NNC*>(c);
	my_continuous_set=nc->get_ptr();
	success = true;
}

void continuous_set_to_PPL_adaptor::dispatch(
		const continuous::support_function::sfm_cont_set<global_types::float_type>* c) {
	my_continuous_set = continuous::continuous_set::ptr();
	success = false;
}

void continuous_set_to_PPL_adaptor::dispatch(
		const continuous::spacetime_flowpipe<global_types::float_type>* c) {
	my_continuous_set = continuous::continuous_set::ptr();
	success = false;
}

void continuous_set_to_PPL_adaptor::dispatch(const continuous::predicate_continuous_set* c) {
	my_continuous_set=ppl_polyhedron::convert_to_continuous_set_PPL_NNC(c->get_predicate());
//std::cerr << " result " << my_continuous_set;
	success = true;
}

void continuous_set_to_PPL_adaptor::dispatch(
		const continuous::continuous_set_simulation<global_types::float_type>* c){
	my_continuous_set = continuous::continuous_set::ptr();
	success = false;
}

void dynamics_to_PPL_adaptor::dispatch(const continuous::ode_affine_dynamics<global_types::rational_type>* c) {
	my_dynamics = continuous::continuous_dynamics::ptr();
	success = false;
}

void dynamics_to_PPL_adaptor::dispatch(const continuous::ode_affine_dynamics<global_types::float_type>* c) {
	my_dynamics = continuous::continuous_dynamics::ptr();
	success = false;
}

transform_to_PPL_adaptor* transform_to_PPL_adaptor::clone() {
	return new transform_to_PPL_adaptor(*this);
}

automaton_to_PPL_adaptor::automaton_to_PPL_adaptor() {
}

automaton_to_PPL_adaptor::~automaton_to_PPL_adaptor() {
}

void automaton_to_PPL_adaptor::init() {
	adapt_continuous_set_visitor_ptr con_ad = adapt_continuous_set_visitor_ptr(
			new continuous_set_to_PPL_adaptor());
	adapt_dynamics_visitor_ptr dyn_ad = adapt_dynamics_visitor_ptr(new dynamics_to_PPL_adaptor());
	dyn_ad->init(con_ad);
	adapt_transform_visitor_ptr trans_ad = adapt_transform_visitor_ptr(
			new transform_to_PPL_adaptor());
	trans_ad->init(con_ad);
	adapt_automaton_visitor::define(con_ad, dyn_ad, trans_ad, global_types::STD_BOOL, global_types::GMP_RATIONAL);
	reset();
}

automaton_to_PPL_adaptor* automaton_to_PPL_adaptor::clone() const {
	return new automaton_to_PPL_adaptor(*this);
}

void automaton_to_PPL_adaptor::epilogue(hybrid_automaton& h) {
	adapt_automaton_visitor::epilogue(h);
//	using std::cout;
//	cout << "PPL adapted " << h.get_name() << ", "<< get_success() << ", ";
//	cout << my_loc_count << " locs, " << my_trans_count << " trans." << std::endl;
}

}
