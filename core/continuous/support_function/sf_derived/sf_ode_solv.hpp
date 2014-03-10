/*
 * sf_ode_solv.hpp
 *
 *  Created on: Oct 11, 2010
 *      Author: frehse
 */

#ifndef SF_ODE_SOLV_HPP_
#define SF_ODE_SOLV_HPP_

#include "sf_ode_solv.h"

#include "math/vdom/vdom_matrix_operators.h"
#include "math/vdom/vdom_vector_operators.h"
#include "math/vdom/affine_map_utility.h"
#include "math/vdom/positional_vdomain.h"

namespace continuous {
namespace support_function {

template<typename scalar_type> sf_ode_solv<scalar_type>::sf_ode_solv(
		const support_function_provider::const_ptr& X0, const typed_dynamics<
				scalar_type>& dyn, scalar_type delta) :
	sf_unary<scalar_type> (s), my_dyn(dyn), my_delta(delta) {
}

template<typename scalar_type> sf_ode_solv<scalar_type>::sf_ode_solv(
		const support_function_provider::const_ptr& X0, const typed_dynamics<
				scalar_type>& dyn, scalar_type delta, const affine_map& M) :
	sf_unary<scalar_type> (s, M), my_dyn(dyn), my_delta(delta) {
}

template<typename scalar_type> sf_ode_solv<scalar_type>::~sf_ode_solv() {
}

template<typename scalar_type> sf_ode_solv<scalar_type>* sf_ode_solv<
		scalar_type>::clone() const {
	support_function_provider::const_ptr new_root(this->my_set->clone());
	if (this->get_map())
		return new sf_ode_solv<scalar_type> (new_root, my_dyn, my_delta,
				*this->get_map());
	else
		return new sf_ode_solv<scalar_type> (new_root, my_dyn, my_delta);
}

template<typename scalar_type> int sf_ode_solv<scalar_type>::get_memory() const {
	throw std::runtime_error("sf_ode_solv : missing implementation get_memory");
	return 0;
}

template<typename scalar_type> continuous_set_predicate::ptr sf_ode_solv<
		scalar_type>::get_predicate() const {
	throw std::runtime_error(
			"sf_ode_solv : missing implementation get_predicate");
	return continuous_set_predicate::ptr();
}

template<typename scalar_type> void sf_ode_solv<scalar_type>::print(
		std::ostream& os) const {
	throw std::runtime_error("sf_ode_solv : missing implementation print");
}

template<typename scalar_type> const variable_id_set& sf_ode_solv<scalar_type>::get_variable_ids() const {
	return sf_unary<scalar_type>::get_variable_ids();
}
template<typename scalar_type> void sf_ode_solv<scalar_type>::reassign_primedness(
		unsigned int, unsigned int) {
	throw std::runtime_error(
			"sf_ode_solv : missing implementation reassign_primedness");
}
template<typename scalar_type> void sf_ode_solv<scalar_type>::increase_primedness(
		unsigned int) {
	throw std::runtime_error(
			"sf_ode_solv : missing implementation increase_primedness");
}
template<typename scalar_type> void sf_ode_solv<scalar_type>::decrease_primedness(
		unsigned int) {
	throw std::runtime_error(
			"sf_ode_solv : missing implementation decrease_primedness");
}

template<typename scalar_type> void sf_ode_solv<scalar_type>::compute_support(
		const math::vdom_vector<Rational>& l, Rational& max_value,
		math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	throw std::runtime_error(
			"sf_ode_solv : missing implementation compute_support Rational");
}

template<typename scalar_type> void sf_ode_solv<scalar_type>::compute_support(
		const math::vdom_vector<double>& l, double& max_value,
		math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) const {

	// define the system for the ode solver
	typed_dynamics_state_functor<scalar_type> f(my_dyn);
	// define the root function
	// compute the trajectory up to delta
	// find the max value

	throw std::runtime_error(
			"sf_ode_solv : missing implementation compute_support double");
}

}
}

#endif /* SF_ODE_SOLV_HPP_ */
