#ifndef ODE_AFFINE_DYNAMICS_H_
#define ODE_AFFINE_DYNAMICS_H_

#include "core/continuous/continuous_dynamics/typed_dynamics.h"

#include "math/vdom/affine_map.h"
#include "math/vdom/index_to_variable_id_map_provider.h"
#include "core/continuous/support_function_provider.h"

namespace continuous {

/** Class to represent the nondeterministic ODE (differential inclusion)
 *  x'==Ax+b+u, where u \in U.
 *
 *	It is assumed that U is "centered" around zero, i.e., that the norm of U is reduced
 *	by including any offset to b. This is, however, not a strict requirement
 *	(U is not required to contain zero).
 *
 *	The domain of u is that of x, i.e., u(x1) is added to x(x1), u(x2) to x(x2) etc.
 */
template<typename scalar_type>
class ode_affine_dynamics: public typed_dynamics<scalar_type> ,
		public math::affine_map<scalar_type> {
public:
	typedef boost::shared_ptr<ode_affine_dynamics<scalar_type> > ptr;
	typedef boost::shared_ptr<const ode_affine_dynamics<scalar_type> >
			const_ptr;
	typedef typename typed_dynamics<scalar_type>::const_visitor const_visitor;

	typedef typename math::affine_map<scalar_type>::matrix_type matrix_type;
	typedef typename math::affine_map<scalar_type>::vector_type vector_type;

	/** Member types */
	typedef typename typed_dynamics<scalar_type>::state state;
	typedef typename typed_dynamics<scalar_type>::matrix matrix;
	typedef support_function_provider offset_type;

	/** Construct empty dynamics */
	ode_affine_dynamics() {};

	/** Construct dynamics dx/dt=A*x+b from an affine map.
	 */
	ode_affine_dynamics(const math::affine_map<scalar_type>& M, const typename offset_type::const_ptr& U_centered=offset_type::const_ptr()) :
		math::affine_map<scalar_type>(M),my_U(U_centered) {
	}
	;

	/** Virtual destructor for possible derived classes */
	virtual ~ode_affine_dynamics() {
	}
	;

	/** Obtain a predicate representation */
	virtual dynamics_predicate::ptr get_predicate() const {
		throw std::runtime_error("no get_relation for ode_affine_dynamics");
	}
	;

	/** Returns the ids of all variables in the domain or codomain.
	 *  A variable x and its derivative \dot x are considered as being the same variable. */
	virtual variable_id_set get_variable_ids() const {
		return get_unprimed_variables(math::affine_map<scalar_type>::get_variable_ids());
	}
	;

	/** Returns the variables that are not in the codomain. */
	virtual variable_id_set get_unconstrained_variable_ids() const {
		variable_id_set allvars=get_variable_ids();
		variable_id_set derivvars=get_primed_variables(allvars);
		variable_id_set statevars=get_unprimed_variables(derivvars);

		variable_id_set algebvars=allvars;
		set_difference_assign(algebvars,statevars);

		return algebvars;
	}
	;

	/** Returns the set of nondeterministic offsets U
	 */
	const offset_type::const_ptr& get_U() const {
		return my_U;
	}
	;

	/** Define a set of nondeterministic offsets
	 *
	 * The domain of U_centered needs to be the same
	 * as that of x.
	 */
	void set_U(const offset_type::const_ptr& U_centered) {
		my_U = U_centered;
	}
	;

	/** Compute dx/dt=f(x).
	 *
	 * If the dynamics are nondeterministic, pick a "representative" point
	 * (e.g., center). */
	state compute_deriv(const state& x) const {
		return map(x);
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
		return this->get_A();
	}
	;

	/** Output as a stream of characters. */
	virtual void print(std::ostream& os) const {
		math::affine_map<scalar_type>::print(os);
		if (get_U()) {
			os << " with offset " << get_U();
		}
	}
	;

	/** Accept a dispatcher. */
	virtual void accept(const_visitor& d) const {
		d.dispatch(this);
	}
	;

private:
	typename offset_type::const_ptr my_U; // nondeteministic offset U
};

}

#endif /*ODE_AFFINE_DYNAMICS_H_*/
