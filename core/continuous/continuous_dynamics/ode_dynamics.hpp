/*
 * ode_dynamics.cpp
 *
 *  Created on: Aug 31, 2009
 *      Author: frehse
 */

#include "core/continuous/continuous_dynamics/ode_dynamics.h"
#include "core/predicates/valuation_function_tree_nodes.h"
#include "core/predicates/valuation_function_tree_ops.h"
#include "utility/shared_ptr_output.h"
#include "utility/tree_node.h"
#include "core/predicates/node_eval_visitor.h"

/** Forward declarations */
namespace tree {
class node;
typedef boost::shared_ptr<node> node_ptr;
}
namespace valuation_functions {
variable_id_set get_unprimed_variable_ids(tree::node_ptr p);
// from #include "core/predicates/valuation_function_tree_ops.h"
template<typename scalar_type>
tree::node::ptr compute_diff(tree::node::ptr p, const variable& var);
}
std::ostream& operator<<(std::ostream& os, const tree::node& p);
//std::ostream& operator<<(std::ostream& os, const tree::node::ptr& p);

namespace continuous {

template<typename scalar_type>
ode_dynamics<scalar_type>::ode_dynamics() :
		has_jacobian(false), has_hessian(false) {
}

template<typename scalar_type>
ode_dynamics<scalar_type>::ode_dynamics(const positional_vdomain& codom,
		const function_vector& assgn) :
		my_codom(codom), my_assgn(assgn), has_jacobian(false), has_hessian(
				false) {
}

template<typename scalar_type>
ode_dynamics<scalar_type>::~ode_dynamics() {
}

template<typename scalar_type>
void ode_dynamics<scalar_type>::insert(variable_id var, const function_ptr& f) {
	my_assgn.push_back(f);
	// extend the domain
	std::vector<variable> vars = my_codom.get_variable_vector();
	vars.push_back(variable(var));
	my_codom = positional_vdomain(vars);
}

template<typename scalar_type>
dynamics_predicate::ptr ode_dynamics<scalar_type>::get_predicate() const {
	throw std::runtime_error("no get_relation for ode_affine_dynamics");
}

template<typename scalar_type>
variable_id_set ode_dynamics<scalar_type>::get_variable_ids() const {
	variable_id_set vis = my_codom.get_variable_ids();
	// for each ode in the list, add the set of variables on the right hand side.
	for (function_vector::const_iterator it = my_assgn.begin();
			it != my_assgn.end(); ++it) {
		// add the set of variables on the right hand side
		variable_id_set v = valuation_functions::get_unprimed_variable_ids(*it);
		vis.insert(v.begin(), v.end());
	}
	return vis;
}

template<typename scalar_type>
variable_id_set ode_dynamics<scalar_type>::get_unconstrained_variable_ids() const {
	variable_id_set allvars = get_variable_ids();

	// get the variables that have a derivative
	variable_id_set statevars = my_codom.get_variable_ids();
	variable_id_set algebvars = allvars;
	set_difference_assign(algebvars, statevars);

	return algebvars;
}

template<typename scalar_type>
typename ode_dynamics<scalar_type>::const_iterator ode_dynamics<scalar_type>::begin() const {
	return my_assgn.begin();
}

template<typename scalar_type>
typename ode_dynamics<scalar_type>::const_iterator ode_dynamics<scalar_type>::end() const {
	return my_assgn.end();
}

template<typename scalar_type>
void ode_dynamics<scalar_type>::print(std::ostream& os) const {
	for (function_vector::size_type i = 0; i < my_assgn.size(); ++i) {
		if (i > 0)
			os << " & ";
		os << my_codom.get_variable(i) << " = " << my_assgn[i];
	}
}

template<typename scalar_type>
void ode_dynamics<scalar_type>::accept(const_visitor& d) const {
	d.dispatch(this);
}

template<typename scalar_type>
typename ode_dynamics<scalar_type>::state ode_dynamics<scalar_type>::compute_deriv(
		const state& x) const {
	// evaluate every one of the functions and return the result
	// with the domain my_codom

	state v_vector(my_codom);

	for (std::vector<function_ptr>::size_type i = 0; i < my_assgn.size(); ++i) {

		v_vector[i] = compute_deriv(i,x);

	}
	return v_vector;
}

template<typename scalar_type>
scalar_type ode_dynamics<scalar_type>::compute_deriv(size_type i, const state& x) const {
	return valuation_functions::arithmetic_eval_node<scalar_type>(my_assgn[i], x);
}


template<typename scalar_type>
const typename ode_dynamics<scalar_type>::function_vector_vector& ode_dynamics<
		scalar_type>::jacobian() const {
	size_t N = my_assgn.size(); // number of equations
	size_t M = my_codom.size(); // number of variables
	if (!has_jacobian) {
		// create and store the jacobian symbolically
		// reserve space
		ode_dynamics<scalar_type>* nonconst_this = const_cast<ode_dynamics<
				scalar_type>*>(this);
		nonconst_this->my_jacobian = function_vector_vector(N,
				function_vector(N));
		for (std::vector<function_ptr>::size_type i = 0; i < N; ++i) {
			for (std::vector<function_ptr>::size_type j = 0; j < M; ++j) {
				nonconst_this->my_jacobian[i][j] =
						valuation_functions::compute_diff<scalar_type>(
								my_assgn[i], my_codom.get_variable(j));
			}
		}
		LOGGER_OS(DEBUG7, __FUNCTION__) << "Jacobian: " << my_jacobian;
		nonconst_this->has_jacobian = true;
	}
	return my_jacobian;
}
;

template<typename scalar_type>
typename ode_dynamics<scalar_type>::matrix ode_dynamics<scalar_type>::compute_jacobian(
		const state& x) const {
	size_t N = my_assgn.size(); // number of equations
	size_t M = my_codom.size(); // number of variables
	const function_vector_vector& jac = jacobian();

	math::vdom_matrix<scalar_type> J(my_codom, my_codom);
	for (std::vector<function_ptr>::size_type i = 0; i < N; ++i) {
		for (std::vector<function_ptr>::size_type j = 0; j < M; ++j) {
			J(i, j) = valuation_functions::arithmetic_eval_node<scalar_type>(
					jac[i][j], x);
		}
	}
	return J;
}
;

template<typename scalar_type>
const typename ode_dynamics<scalar_type>::function_vector_vector_vector& ode_dynamics<
		scalar_type>::hessian() const {
	size_t N = my_assgn.size(); // number of equations
	size_t M = my_codom.size(); // number of variables
	if (!has_hessian) {
		// create and store the jacobian symbolically
		// note: because of symmetry, only upper triangle is computed
		// reserve space
		ode_dynamics<scalar_type>* nonconst_this = const_cast<ode_dynamics<
				scalar_type>*>(this);
		nonconst_this->my_hessian = function_vector_vector_vector(N,
				function_vector_vector(N, function_vector(N)));
		for (std::vector<function_ptr>::size_type i = 0; i < N; ++i) {
			for (std::vector<function_ptr>::size_type j = 0; j < M; ++j) {
				function_ptr dfi_dxj = jacobian()[i][j];
						// valuation_functions::compute_diff<scalar_type>(my_assgn[i], my_codom.get_variable(j));

				for (std::vector<function_ptr>::size_type k = j; k < M; ++k) {
					// debug output
					function_ptr dfi_dxj_dxk =
							valuation_functions::compute_diff<scalar_type>(
									dfi_dxj, my_codom.get_variable(k));
					nonconst_this->my_hessian[i][j][k] = dfi_dxj_dxk;
				}
			}
		}
		LOGGER_OS(DEBUG7, __FUNCTION__) << "Hessian: " << my_hessian;

		nonconst_this->has_hessian = true;
	}
	return my_hessian;
}
;

template<typename scalar_type>
typename std::vector<typename continuous::ode_dynamics<scalar_type>::matrix> ode_dynamics<
		scalar_type>::compute_hessian(const state& x) const {
	size_t N = my_assgn.size(); // number of equations
	size_t M = my_codom.size(); // number of variables
	const function_vector_vector_vector& hess = hessian();

	std::vector<math::vdom_matrix<scalar_type> > J(N,
			math::vdom_matrix<scalar_type>(my_codom, my_codom));
	for (std::vector<function_ptr>::size_type i = 0; i < N; ++i) {
		for (std::vector<function_ptr>::size_type j = 0; j < M; ++j) {
			for (std::vector<function_ptr>::size_type k = j; k < M; ++k) {
				J[i](j, k) = valuation_functions::arithmetic_eval_node<
						scalar_type>(hess[i][j][k], x);
				J[i](k, j) = J[i](j, k); // exploit symmetry
			}
		}
	}
	return J;
}
;

template<typename scalar_type>
scalar_type ode_dynamics<scalar_type>::compute_hessian(size_type i, size_type j, size_type k, const state& x) const {
	const function_vector_vector_vector& hess = hessian();
	if (k<j) {
		std::swap(k,j);
	}
	scalar_type val = valuation_functions::arithmetic_eval_node<scalar_type>(hess[i][j][k], x);
	return val;
}
;

template<typename scalar_type>
const positional_vdomain& ode_dynamics<scalar_type>::codom() const {
	return my_codom;
}

}

