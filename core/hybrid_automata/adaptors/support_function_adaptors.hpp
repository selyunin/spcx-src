/*
 * support_function_adapters.hpp
 *
 *  Created on: Nov 9, 2009
 *      Author: frehse
 */

#ifndef SUPPORT_FUNCTION_ADAPTERS_HPP_
#define SUPPORT_FUNCTION_ADAPTERS_HPP_

//#include "support_function_adaptors.h"
#include "core/continuous/continuous_dynamics/continuous_dynamics.h"
#include "core/continuous/continuous_set_transforms/continuous_set_transforms.h"
#include "core/continuous/predicate_continuous_set.h"
#include "core/continuous/convert_predicate.h"
#include "core/continuous/support_function/sf_base/sf_unary.h"
#include "core/continuous/polyhedra/hyperbox/bounding_box.h"

#include "core/hybrid_automata/hybrid_automaton.h"

#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_operators.h"
#include "core/continuous/polyhedra/constr_polyhedron/constr_polyhedron_constructors.h"
#include "core/hybrid_automata/assignment_traitment.h"

#include "core/continuous/support_function_provider.h"

#include "math/ode_solving/traj_simu/continuous_set_simulation.h"

#include "math/vdom/affine_map.h"

namespace hybrid_automata {

// -------------------------------------------------------------------
// Implementations
// -------------------------------------------------------------------

template<typename scalar_type>
void continuous_set_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::support_function_provider* c) {
	my_continuous_set = continuous::continuous_set::ptr();
	success = false;
}

template<typename scalar_type>
void continuous_set_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::constr_polyhedron<Rational>* c) {
	my_continuous_set = continuous::construct_constr_polyhedron<scalar_type,
			Rational>(*c);
	success = true;
}

template<typename scalar_type>
void continuous_set_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::constr_polyhedron<double>* c) {
	my_continuous_set = continuous::construct_constr_polyhedron<scalar_type,
			double>(*c);
	success = true;
}

template<typename scalar_type>
void continuous_set_to_supp_f_adaptor<scalar_type>::dispatch(
		const ppl_polyhedron::continuous_set_PPL_NNC* c) {
	my_continuous_set = continuous::construct_constr_polyhedron<scalar_type,
			Rational>((const continuous::polyhedron<Rational>&) (*c));
	success = true;
}

template<typename scalar_type>
void continuous_set_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::support_function::sfm_cont_set<
				global_types::float_type>* c) {
	my_continuous_set = continuous::continuous_set::ptr();
	success = false;
}

template<typename scalar_type>
void continuous_set_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::spacetime_flowpipe<
				global_types::float_type>* c) {
	my_continuous_set = continuous::continuous_set::ptr();
	success = false;
}

template<typename scalar_type>
void continuous_set_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::predicate_continuous_set* c) {
	tree::node::ptr p = c->get_predicate();
	my_continuous_set = continuous::construct_constr_polyhedron<scalar_type>(p);
	success = true;
}

template<typename scalar_type>
void continuous_set_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::continuous_set_simulation<global_types::float_type>* c){
	my_continuous_set = continuous::continuous_set::ptr();
	success = false;
}

template<typename scalar_type>
void dynamics_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::relation_dynamics* d) {
	using namespace continuous;

	continuous::relation::const_ptr c = d->get_relation();

	// convert the relation using continuous_set_adaptor
	tree::node::ptr p = c->get_predicate();

	continuous::continuous_dynamics::ptr dyn;
	math::affine_map<scalar_type> M;

	//std::cout << "adapting with my_accept_nonlinear:" << std::boolalpha << my_accept_nonlinear << std::endl;
	try {
		// convert to affine_map, requiring that the primed variables are also part of the domain
		success = continuous::convert_predicate_to_affine_map<scalar_type>(p,
				M, true);
	} catch (std::exception& e) {
		throw basic_exception("The following are not affine dynamics:\n"
				+ logger::formatted_to_string(p), e);
	}
	//	std::cout << M.get_A() << M.get_b() << std::endl;
	//	std::cout << M.get_A().domain() << M.get_A().codomain() << M.get_b().domain() << std::endl;
	//	std::cout << M.get_A().get_matrix() << M.get_b().get_vector() << std::endl;
	if (success) {
		//		std::cout << "converted dynamics " << p << " to " << M << std::endl;
		dyn = typename continuous::ode_affine_dynamics<scalar_type>::ptr(
				new continuous::ode_affine_dynamics<scalar_type>(M));
	} else if (my_accept_nonlinear) {
		positional_vdomain codom; // we will obtain the codomain of the function vector
		std::vector<tree::node::ptr> funcs; // the vector of functions
		success = convert_predicate_to_function_vector<scalar_type>(p, codom, funcs);
		if (success) {
			// std::cout << "converted dynamics " << p << " to " << M << std::endl;
			dyn = typename continuous::ode_dynamics<scalar_type>::ptr(
					new continuous::ode_dynamics<scalar_type>(codom,funcs));
		} else {
			throw basic_exception("The following are not ode dynamics:\n"
					+ logger::formatted_to_string(p));
		}
	} else {
		throw basic_exception("The following are not affine dynamics:\n"
				+ logger::formatted_to_string(p));
	}

	my_dynamics = dyn;
}

template<typename scalar_type>
void dynamics_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::ode_affine_dynamics<global_types::rational_type>* c) {
	my_dynamics = typename continuous::ode_affine_dynamics<scalar_type>::ptr(
			new continuous::ode_affine_dynamics<scalar_type>(c->convert_to<
					scalar_type> ()));
	success = true;
}

template<typename scalar_type>
void dynamics_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::ode_affine_dynamics<global_types::float_type>* c) {
	my_dynamics = typename continuous::ode_affine_dynamics<scalar_type>::ptr(
			new continuous::ode_affine_dynamics<scalar_type>(c->convert_to<
					scalar_type> ()));
	success = true;
}

template<typename scalar_type>
transform_to_supp_f_adaptor<scalar_type>* transform_to_supp_f_adaptor<
		scalar_type>::clone() {
	return new transform_to_supp_f_adaptor(*this);
}

template<typename scalar_type>
void transform_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::relation_transform* t) {
	continuous::relation::const_ptr c = t->get_relation_const();
	tree::node::ptr p = c->get_predicate();

	// convert the relation using continuous_set_adaptor

	bool aff_dyn = false;
	continuous::continuous_set_transform::ptr transf;
	math::affine_map<scalar_type> M;

	try {
		success
				= continuous::convert_predicate_to_affine_map<scalar_type>(p, M);
	} catch (std::exception& e) {
		std::stringstream s;
		logger::copyfmt_to(s);
		s << p;
		throw basic_exception("The following is not an affine transform:\n"
				+ s.str(), e);
	}
	//std::cout << "converted transform " << p << " to " << M << std::endl;
	//std::cout << M.get_A() << M.get_b() << std::endl;
	//std::cout << M.get_A().domain() << M.get_A().codomain() << M.get_b().domain() << std::endl;
	//std::cout << M.get_A().get_matrix() << M.get_b().get_vector() << std::endl;
	if (success) {
		transf = continuous::continuous_set_transform::ptr(
				new continuous::reset_affine_transform<scalar_type>(M));
	} else {
		std::stringstream s;
		logger::copyfmt_to(s);
		s << p;
		throw basic_exception("The following is not an affine transform:\n"
				+ s.str());
	}
	my_transform = transf;
}

template<typename scalar_type>
void transform_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::reset_affine_transform<global_types::rational_type>* t) {
	my_transform = continuous::continuous_set_transform::ptr(
			new continuous::reset_affine_transform<scalar_type>(t->convert_to<
					scalar_type> ()));
	success = true;
}

template<typename scalar_type>
void transform_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::reset_affine_transform<global_types::float_type>* t) {
	my_transform = continuous::continuous_set_transform::ptr(
			new continuous::reset_affine_transform<scalar_type>(t->convert_to<
					scalar_type> ()));
	success = true;
}

template<typename scalar_type>
void transform_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::reset_function_transform* t) {
	// @todo to be done
	success = false;
}

template<typename scalar_type>
void transform_to_supp_f_adaptor<scalar_type>::dispatch(
		const continuous::constant_bound_time_elapse_transform* t) {
	// @todo to be done
	success = false;
}

inline automaton_to_supp_f_adaptor::automaton_to_supp_f_adaptor(
		coefficient_type bool_t, coefficient_type number_t, bool accept_nonlinear) {
	set_bool_type(bool_t);
	set_number_type(number_t);
	my_accept_nonlinear = accept_nonlinear;
}

inline automaton_to_supp_f_adaptor::~automaton_to_supp_f_adaptor() {
}

inline void automaton_to_supp_f_adaptor::init() {
	adapt_continuous_set_visitor_ptr con_ad = adapt_continuous_set_visitor_ptr(
			global_types::template_type_factory<adapt_continuous_set_visitor,
					continuous_set_to_supp_f_adaptor>::create_number(
					get_number_type()));
	adapt_dynamics_visitor_ptr dyn_ad = adapt_dynamics_visitor_ptr(
			global_types::template_type_factory<adapt_dynamics_visitor,
					dynamics_to_supp_f_adaptor>::create_number(
					get_number_type(),my_accept_nonlinear));
	dyn_ad->init(con_ad);
	adapt_transform_visitor_ptr trans_ad = adapt_transform_visitor_ptr(
			global_types::template_type_factory<adapt_transform_visitor,
					transform_to_supp_f_adaptor>::create_number(
					get_number_type()));
	trans_ad->init(con_ad);
	adapt_automaton_visitor::define(con_ad, dyn_ad, trans_ad, get_bool_type(),
			get_number_type());
}

inline void automaton_to_supp_f_adaptor::epilogue(hybrid_automaton& h) {
	adapt_automaton_visitor::epilogue(h);
	//	using std::cout;
	//	cout << "PPL adapted " << h.get_name() << ", "<< get_success() << ", ";
	//	cout << my_loc_count << " locs, " << my_trans_count << " trans." << std::endl;
}

}

#endif /* SUPPORT_FUNCTION_ADAPTERS_HPP_ */
