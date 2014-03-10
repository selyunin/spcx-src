/*
 * typed_dynamics.h
 *
 *  Created on: Oct 7, 2010
 *      Author: frehse
 */

#ifndef TYPED_DYNAMICS_H_
#define TYPED_DYNAMICS_H_

#include "core/continuous/continuous_dynamics/continuous_dynamics_base.h"
#include "math/ode_solving/ode_solver.h"

namespace continuous {

/** Abstract interface class for representing set-based continuous dynamics, i.e.,
 * set based ODEs of the form \f$ \dot x \in F(x) \f$.
 *
 * This derivation of continuous_dynamics offers the function f as deriv(...).
 */
template<typename scalar_type>
class typed_dynamics: public continuous_dynamics {
public:
	typedef boost::shared_ptr<typed_dynamics<scalar_type> > ptr;
	typedef boost::shared_ptr<const typed_dynamics<scalar_type> > const_ptr;

	typedef typename math::vdom_vector<scalar_type>::size_type size_type;
	typedef math::vdom_vector<scalar_type> state;
	typedef math::vdom_matrix<scalar_type> matrix;

	/** Return a shared_ptr to *this. */
	ptr get_ptr() {
		return continuous_dynamics::get_ptr();
	}
	;

	/** Return a shared_ptr to const *this. */
	const_ptr get_const_ptr() const {
		return continuous_dynamics::get_const_ptr();
	}
	;

	virtual ~typed_dynamics() {
	}
	;

	/** Compute the derivative of x_i at state x.
	 *
	 * The derivative is defined by dx_i/dt=f(x).
	 *
	 * If the dynamics are nondeterministic, pick a "representative" point
	 * (e.g., center). */
	virtual scalar_type compute_deriv(size_type i, const state& x) const {
		throw std::runtime_error(
				"compute_deriv not implemented for this type of dynamics");
		return scalar_type();
	}
	;

	/** Compute the derivative at state x.
	 *
	 * The derivative is defined by dx/dt=f(x).
	 *
	 * If the dynamics are nondeterministic, pick a "representative" point
	 * (e.g., center). */
	virtual state compute_deriv(const state& x) const {
		throw std::runtime_error(
				"compute_deriv not implemented for this type of dynamics");
		return state();
	}
	;

	/** Compute the Jacobian at state x.
	 *
	 * The Jacobian is a matrix where the element (i,j) is defined by
	 * df_i(x)/dx_j.
	 *
	 * If the dynamics are nondeterministic, pick a "representative" point
	 * (e.g., center). */
	virtual matrix compute_jacobian(const state& x) const {
		throw std::runtime_error(
				"compute_jacobian not implemented for this type of dynamics");
		return matrix();
	}
	;

	/** Compute the Hessian of variable x_i at state x with respect to x_j and x_k.
	 *
	 * The Hessian of d/dt x_i = f_i(x) is defined by
	 * d^2 f_i(x)/(dx_j,dx_k).
	 *
	 * If the dynamics are nondeterministic, pick a "representative" point
	 * (e.g., center). */
	virtual scalar_type compute_hessian(size_type i, size_type j, size_type k, const state& x) const {
		throw std::runtime_error(
				"compute_hessian not implemented for this type of dynamics");
		return scalar_type();
	}
	;

	/** Returns the codomain */
//#pragma GCC diagnostic push
//#pragma GCC diagnostic ignored "-Wreturn-type"
	virtual const positional_vdomain& codom() const {
		throw std::runtime_error(
				"positional_vdomain not implemented for this type of dynamics");
//		return positional_vdomain();
	}
	// turn the warnings back on
//#pragma GCC diagnostic pop

};

/** An interface class that allows to use any typed dynamics object as a state functor */
template<typename scalar_type>
class typed_dynamics_state_functor: public math::state_functor<scalar_type> {
public:
	/** Member types */
	typedef typename math::state_functor<scalar_type>::state state;

	/** Maps x to f(x).
	 *
	 * @note We use shared pointers because a reference could end up dangling
	 * if the passed dynamics object is deleted. */
	typed_dynamics_state_functor(
			typename typed_dynamics<scalar_type>::const_ptr dyn) :
		my_dyn(dyn) {
	}
	;
	virtual ~typed_dynamics_state_functor() {
	}
	;
	virtual state map(const state& x) const {
		return my_dyn->compute_deriv(x);
	}
	;
private:
	typename typed_dynamics<scalar_type>::const_ptr my_dyn;
};

/** An interface class that returns a state_matrix_functor that computes the
 * Jacobian for the dynamics */
template<typename scalar_type>
class typed_dynamics_jacobian_functor: public math::state_matrix_functor<scalar_type> {
public:
	/** Member types */
	typedef typename math::state_matrix_functor<scalar_type>::state state;
	typedef typename math::state_matrix_functor<scalar_type>::matrix matrix;

	/** Maps x to f(x).
	 *
	 * @note We use shared pointers because a reference could end up dangling
	 * if the passed dynamics object is deleted. */
	typed_dynamics_jacobian_functor(
			typename typed_dynamics<scalar_type>::const_ptr dyn) :
		my_dyn(dyn) {
	}
	;
	virtual ~typed_dynamics_jacobian_functor() {
	}
	;
	virtual matrix map(const state& x) const {
		return my_dyn->compute_jacobian(x);
	}
	;
private:
	typename typed_dynamics<scalar_type>::const_ptr my_dyn;
};

}

#endif /* TYPED_DYNAMICS_H_ */
