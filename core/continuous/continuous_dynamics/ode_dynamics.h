#ifndef ODE_DYNAMICS_H_
#define ODE_DYNAMICS_H_

#include "core/continuous/continuous_dynamics/continuous_dynamics_base.h"

namespace continuous {

/** Deterministic ODE (Ordinary Differential Equation) system, i.e.,
 *  \f$ \dot x_i = f_i(x_i) \f$ over some \f$i\f$. */
template<typename scalar_type>
class ode_dynamics : public typed_dynamics<scalar_type> {
public:
	typedef boost::shared_ptr<ode_dynamics<scalar_type> > ptr;
	typedef boost::shared_ptr<const ode_dynamics<scalar_type> >
			const_ptr;
	typedef typename typed_dynamics<scalar_type>::const_visitor const_visitor;

	typedef tree::node::ptr function_ptr;
	typedef std::vector<function_ptr> function_vector;
	typedef function_vector::const_iterator const_iterator;
	typedef std::vector<std::vector<function_ptr> > function_vector_vector;
	typedef std::vector<std::vector<std::vector<function_ptr> > > function_vector_vector_vector;

	/** Member types */
	typedef typename typed_dynamics<scalar_type>::size_type size_type;
	typedef typename typed_dynamics<scalar_type>::state state;
	typedef typename typed_dynamics<scalar_type>::matrix matrix;

	/** Construct empty ODE dynamics
	 */
	ode_dynamics();

	/** Construct ODE dynamics from a given codomain a vector of functions
	 *
	 * The codomain is the domain of the derivatives \dot x_i to which
	 * the function funcs[i] is attributed.
	 *
	 * Note that x_i are the unprimed variables.
	 */
	ode_dynamics(const positional_vdomain& codom, const function_vector& funcs);

	/** Virtual destructor for possible derived classes */
	virtual ~ode_dynamics();

	/** Add an ODE.
	 *
	 * Adds the ODE \dot var = f.
	 * Note that var is the unprimed variable.
	 */
	void insert(variable_id var, const function_ptr& f);

	/** Obtain a predicate representation */
	virtual dynamics_predicate::ptr get_predicate() const ;

	/** Returns the ids of all variables in the domain or codomain.
	 *  A variable x and its derivative \dot x are considered as being the same variable. */
	variable_id_set get_variable_ids() const;

	/** Returns the variables whose derivative is not assigned. */
	virtual variable_id_set get_unconstrained_variable_ids() const;

	/** Compute the derivative of x_i at state x.
	 *
	 * The derivative is defined by dx_i/dt=f(x).
	 *
	 * If the dynamics are nondeterministic, pick a "representative" point
	 * (e.g., center). */
	virtual scalar_type compute_deriv(size_type i, const state& x) const;

	/** Compute dx/dt=f(x).
	 *
	 * If the dynamics are nondeterministic, pick a "representative" point
	 * (e.g., center). */
	state compute_deriv(const state& x) const;

	/** Compute the Jacobian at state x.
	 *
	 * The Jacobian is a matrix where the element (i,j) is defined by
	 * df_i(x)/dx_j.
	 *
	 * If the dynamics are nondeterministic, pick a "representative" point
	 * (e.g., center). */
	matrix compute_jacobian(const state& x) const;

	/** Compute the Hessian at state x.
	 *
	 * The Hessian is a vector of matrices where the element (i,j,k) is defined by
	 * d^2 f_i(x)/(dx_j,dx_k).
	 *
	 * If the dynamics are nondeterministic, pick a "representative" point
	 * (e.g., center). */
	std::vector<matrix> compute_hessian(const state& x) const;

	/** Compute the Jacobian as vector-vector of functions
	 *
	 * The Jacobian is a matrix where the element (i,j) is defined by
	 * df_i(x)/dx_j.
	 *
	 */
	const function_vector_vector& jacobian() const;

	/** Compute the Hessian as a vector-vector-vector of functions.
	 *
	 * The Hessian is a vector of matrices where the element (i,j,k) is defined by
	 * d^2 f_i(x)/(dx_j,dx_k).
	 *
	 */
	const function_vector_vector_vector& hessian() const;

	/** Compute the Hessian of variable x_i at state x with respect to x_j and x_k.
	 *
	 * The Hessian of d/dt x_i = f_i(x) is defined by
	 * d^2 f_i(x)/(dx_j,dx_k).
	 *
	 * If the dynamics are nondeterministic, pick a "representative" point
	 * (e.g., center). */
	scalar_type compute_hessian(size_type i, size_type j, size_type k, const state& x) const;

	/** Returns the codomain */
	const positional_vdomain& codom() const;

	const_iterator begin() const;
	const_iterator end() const;
	virtual void print(std::ostream& os) const;

	/** Accept a dispatcher. */
	virtual void accept(const_visitor& d) const;

private:
	positional_vdomain my_codom;
	function_vector my_assgn;
	bool has_jacobian;
	bool has_hessian;
	function_vector_vector my_jacobian;
	function_vector_vector_vector my_hessian;
};

}

#include "ode_dynamics.hpp"

#endif /*ODE_DYNAMICS_H_*/
