/*
 * constr_poly_adaptors.cpp
 *
 *  Created on: Sep 11, 2009
 *      Author: frehse
 */

#include "core/hybrid_automata/adaptors/constr_poly_adaptors.h"

#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_constructors.h"
#include "core/continuous/polyhedra/ppl_polyhedron/continuous_set_PPL_NNC.h"

#include "core/continuous/continuous_dynamics/continuous_dynamics.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transforms.h"
#include "core/continuous/predicate_continuous_set.h"

#include "math/ode_solving/traj_simu/continuous_set_simulation.h"

#include "core/hybrid_automata/hybrid_automaton.h"

namespace hybrid_automata {

continuous_set_to_constr_poly_adaptor::continuous_set_to_constr_poly_adaptor(
		coefficient_type bool_t, coefficient_type number_t) :
	my_bool_type(bool_t), my_number_type(number_t) {
}

void continuous_set_to_constr_poly_adaptor::dispatch(const continuous::support_function_provider* c) {
	my_continuous_set = continuous::continuous_set_ptr();
	success = false;
}

void continuous_set_to_constr_poly_adaptor::dispatch(
		const continuous::constr_polyhedron<Rational>* c) {
	if (my_number_type == global_types::GMP_RATIONAL) {
		// Nothing to be done
		// @todo this const-cast is dirty, but it saves us from copying the object
		continuous::constr_polyhedron<Rational>* nc = const_cast<continuous::constr_polyhedron<
				Rational>*> (c);
		my_continuous_set = nc->get_ptr();
	} else {
		my_continuous_set = continuous::construct_constr_polyhedron(*c, my_number_type);
	}
	success = true;
}
void continuous_set_to_constr_poly_adaptor::dispatch(const continuous::constr_polyhedron<double>* c) {
	if (my_number_type == global_types::STD_DOUBLE) {
		// Nothing to be done
		// @todo this const-cast is dirty, but it saves us from copying the object
		continuous::constr_polyhedron<double>* nc =
				const_cast<continuous::constr_polyhedron<double>*> (c);
		my_continuous_set = nc->get_ptr();
	} else {
		my_continuous_set = continuous::construct_constr_polyhedron(*c, my_number_type);
	}
	success = true;
}
void continuous_set_to_constr_poly_adaptor::dispatch(
		const ppl_polyhedron::continuous_set_PPL_NNC* c) {
	my_continuous_set = continuous::construct_constr_polyhedron(*c, my_number_type);
	success = true;
}

void continuous_set_to_constr_poly_adaptor::dispatch(const continuous::support_function::sfm_cont_set<global_types::float_type>* c) {
	my_continuous_set = continuous::continuous_set_ptr();
	success = false;
}

void continuous_set_to_constr_poly_adaptor::dispatch(const continuous::spacetime_flowpipe<global_types::float_type>* c) {
	my_continuous_set = continuous::continuous_set_ptr();
	success = false;
}

void continuous_set_to_constr_poly_adaptor::dispatch(const continuous::predicate_continuous_set* c) {
	my_continuous_set = continuous::construct_constr_polyhedron(c->get_predicate(), my_number_type);
	success = true;
}

void continuous_set_to_constr_poly_adaptor::dispatch(
		const continuous::continuous_set_simulation<global_types::float_type>* c){
	my_continuous_set = continuous::continuous_set::ptr();
	success = false;
}

void dynamics_to_constr_poly_adaptor::dispatch(const continuous::ode_affine_dynamics<global_types::rational_type>* c) {
	success = false;
}

void dynamics_to_constr_poly_adaptor::dispatch(const continuous::ode_affine_dynamics<global_types::float_type>* c) {
	success = false;
}

transform_to_constr_poly_adaptor* transform_to_constr_poly_adaptor::clone() {
	return new transform_to_constr_poly_adaptor(*this);
}

automaton_to_constr_poly_adaptor::automaton_to_constr_poly_adaptor(coefficient_type bool_t,
		coefficient_type number_t) {
	set_bool_type(bool_t);
	set_number_type(number_t);
}

automaton_to_constr_poly_adaptor::~automaton_to_constr_poly_adaptor() {
}

void automaton_to_constr_poly_adaptor::init() {
	adapt_continuous_set_visitor_ptr con_ad = adapt_continuous_set_visitor_ptr(
			new continuous_set_to_constr_poly_adaptor(get_bool_type(), get_number_type()));
	adapt_dynamics_visitor_ptr dyn_ad = adapt_dynamics_visitor_ptr(
			new dynamics_to_constr_poly_adaptor());
	dyn_ad->init(con_ad);
	adapt_transform_visitor_ptr trans_ad = adapt_transform_visitor_ptr(
			new transform_to_constr_poly_adaptor());
	trans_ad->init(con_ad);
	adapt_automaton_visitor::define(con_ad, dyn_ad, trans_ad, get_bool_type(), get_number_type());
	reset();
}

automaton_to_constr_poly_adaptor* automaton_to_constr_poly_adaptor::clone() const {
	return new automaton_to_constr_poly_adaptor(*this);
}

void automaton_to_constr_poly_adaptor::epilogue(hybrid_automaton& h) {
	adapt_automaton_visitor::epilogue(h);
	//	using std::cout;
	//	cout << "constr_poly adapted " << h.get_name() << ", "<< get_success() << ", ";
	//	cout << my_loc_count << " locs, " << my_trans_count << " trans." << std::endl;
}

}
