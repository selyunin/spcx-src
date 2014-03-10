/*
 * sf_ode_deterministic.hpp
 *
 *  Created on: Oct 11, 2010
 *      Author: frehse
 */

#ifndef sf_ode_deterministic_HPP_
#define sf_ode_deterministic_HPP_

#include "sf_ode_deterministic.h"

#include "math/vdom/vdom_matrix_operators.h"
#include "math/vdom/vdom_vector_operators.h"
#include "math/vdom/affine_map_utility.h"
#include "math/vdom/positional_vdomain.h"
#include "math/ode_solving/ode_solver_chooser.h"
#include "math/vdom/state_vector_functor_utility.h"

namespace continuous {
namespace support_function {

template<typename scalar_type> sf_ode_deterministic<scalar_type>::sf_ode_deterministic(
		const state& x0, typename typed_dynamics<scalar_type>::const_ptr dyn,
		scalar_type delta) :
	my_x0(x0), my_dyn(dyn), my_delta(delta) {
	my_solver = math::ode::ode_solver_chooser<scalar_type>::get_instance();
}

template<typename scalar_type> sf_ode_deterministic<scalar_type>::sf_ode_deterministic(
		const state& x0, typename typed_dynamics<scalar_type>::const_ptr dyn,
		scalar_type delta, const affine_map& M) :
	sf_set<scalar_type> (M), my_x0(x0), my_dyn(dyn), my_delta(delta) {
	my_solver = math::ode::ode_solver_chooser<scalar_type>::get_instance();
}

template<typename scalar_type> sf_ode_deterministic<scalar_type>::~sf_ode_deterministic() {
	delete my_solver;
}

template<typename scalar_type> sf_ode_deterministic<scalar_type>* sf_ode_deterministic<
		scalar_type>::clone() const {
	if (this->get_map())
		return new sf_ode_deterministic<scalar_type> (my_x0, my_dyn, my_delta,
				*this->get_map());
	else
		return new sf_ode_deterministic<scalar_type> (my_x0, my_dyn, my_delta);
}

template<typename scalar_type>
math::tribool sf_ode_deterministic<scalar_type>::is_universe() const {
	return false;
}

template<typename scalar_type>
math::tribool sf_ode_deterministic<scalar_type>::is_empty() const {
	return false;
}

template<typename scalar_type> int sf_ode_deterministic<scalar_type>::get_memory() const {
	throw std::runtime_error(
			"sf_ode_deterministic : missing implementation get_memory");
	return 0;
}

template<typename scalar_type> continuous_set_predicate::ptr sf_ode_deterministic<
		scalar_type>::get_predicate() const {
	throw std::runtime_error(
			"sf_ode_deterministic : missing implementation get_predicate");
	return continuous_set_predicate::ptr();
}

template<typename scalar_type> void sf_ode_deterministic<scalar_type>::print(
		std::ostream& os) const {
	throw std::runtime_error(
			"sf_ode_deterministic : missing implementation print");
}

template<typename scalar_type> const variable_id_set& sf_ode_deterministic<
		scalar_type>::get_variable_ids() const {
	if (!this->get_map() || this->get_map()->is_empty())
		return my_x0.get_variable_ids();
	else
		return this->get_map()->codomain().get_variable_ids();
}
template<typename scalar_type> void sf_ode_deterministic<scalar_type>::reassign_primedness(
		unsigned int, unsigned int) {
	throw std::runtime_error(
			"sf_ode_deterministic : missing implementation reassign_primedness");
}
template<typename scalar_type> void sf_ode_deterministic<scalar_type>::increase_primedness(
		unsigned int) {
	throw std::runtime_error(
			"sf_ode_deterministic : missing implementation increase_primedness");
}
template<typename scalar_type> void sf_ode_deterministic<scalar_type>::decrease_primedness(
		unsigned int) {
	throw std::runtime_error(
			"sf_ode_deterministic : missing implementation decrease_primedness");
}

template<typename scalar_type> void sf_ode_deterministic<scalar_type>::compute_support(
		const math::vdom_vector<Rational>& l, Rational& max_value,
		math::vdom_vector<Rational>& support_vec, bool& is_empty,
		bool& is_bounded) const {
	throw std::runtime_error(
			"sf_ode_deterministic : missing implementation compute_support Rational");
}

template<typename scalar_type> void sf_ode_deterministic<scalar_type>::compute_support(
		const math::vdom_vector<double>& l, double& max_value,
		math::vdom_vector<double>& support_vec, bool& is_empty,
		bool& is_bounded) const {

	// define the system for the ode solver
	typed_dynamics_state_functor<scalar_type> f(my_dyn);
	// define the root function
	// note: we need to remember to add the inhomogenuous term of l to the result
	//       for now we compute the max without the inh. term (using lvec)
	math::vdom_vector<scalar_type> lvec = l.template convert_to<scalar_type> ();
	//typed_dynamics_scalar_product_functor<scalar_type> g(my_dyn, lvec);

	// construct the vector of arguments for the scalar products
	typedef typename math::scalar_product_functor<scalar_type>::state state;
	typename math::scalar_product_functor<scalar_type>::state_vector V(1);
	V[0] = lvec;
	math::scalar_product_functor<scalar_type> g(f, V);

	// compute the trajectory up to delta
	my_solver->init_rootf(my_x0, scalar_type(0), f, g);
	typename math::ode::ode_solver<scalar_type>::rootf_result res;
	res = my_solver->solve_rootf(my_delta);

	// find the max value:
	scalar_type v;
	// it's either a) x0, b) res.stop_state, or c) one of the roots
	// a)
	math::vdom_vector<scalar_type> sfvec(my_x0);
	scalar_type sfval = scalar_product(my_x0, lvec);
	// b)
	v = scalar_product(res.stop_state, lvec);
	if (v > sfval) {
		sfval = v;
		sfvec = res.stop_state;
	}
	// c)
	if (res.traj.size() > 0) {
//std::cout << "roots at states:" << std::endl << res.traj << std::endl;
		math::vector<scalar_type> roots_sfval = res.traj.get_states()
				* lvec.get_vector();
		unsigned int max_index;
		v = max(roots_sfval, max_index);
		if (v > sfval) {
			sfval = v;
			sfvec = res.traj.get_state(max_index);
		}
	}

	is_empty = false;
	is_bounded = true;
	max_value = convert_element<double> (sfval);
	support_vec = sfvec.template convert_to<double> ();
}

}
}

#endif /* sf_ode_deterministic_HPP_ */
