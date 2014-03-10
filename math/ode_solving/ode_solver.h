/*
 * ode_solver.h
 *
 *  Created on: Oct 7, 2010
 *      Author: frehse
 */

#ifndef ODE_SOLVER_H_
#define ODE_SOLVER_H_

#include <boost/shared_ptr.hpp>
#include "math/vdom/affine_map.h"
#include "math/vdom/trajectory.h"
#include "math/vdom/state_functor.h"
#include "math/vdom/state_vector_functor.h"
#include "math/vdom/state_matrix_functor.h"

#include "math/vdom/lin_constraint_system.h"
#include <sstream>

namespace math {
namespace ode {

/** A class for storing the default parameters *
 */
class ode_defaults {
public:
	static double get_abs_tol();
	static double get_rel_tol();
	static void set_abs_tol(double x);
	static void set_rel_tol(double x);
private:
	static double my_abs_tol;
	static double my_rel_tol;
};

/** An interface for ode solvers.
 * scalar_type must support the instantiations scalar_type(0) and scalar_type(1).
 */
template<typename scalar_type> class ode_solver {
public:
	typedef boost::shared_ptr<ode_solver<scalar_type> > ptr;
	typedef boost::shared_ptr<const ode_solver<scalar_type> > const_ptr;

	/** Classes for defining dynamics */
	typedef math::state_functor<scalar_type> state_functor;
	typedef typename state_functor::state state;

	/** Classes for defining root functions */
	typedef math::state_vector_functor<scalar_type> state_vector_functor;
	typedef typename state_vector_functor::vector vector;

	/** Classes for defining Jacobian informations */
	typedef math::state_matrix_functor<scalar_type> state_matrix_functor;
	typedef typename state_matrix_functor::matrix matrix;

	/** Virtual Destructor. */
	virtual ~ode_solver();

	struct ode_parameters {
		/** Default parameters. */
		ode_parameters() :
			max_timestep(-1.0), rel_tol(ode_defaults::get_rel_tol()), abs_tol(
					ode_defaults::get_abs_tol()) {
		}
		;
		template<typename result_type>
		typename ode_solver<result_type>::ode_parameters convert_to() {
			typename ode_solver<result_type>::ode_parameters par;
			par.max_timestep = max_timestep;
			par.rel_tol = rel_tol;
			par.abs_tol = abs_tol;
			return par;
		}
		;
		/** Max time step of the solver (negative for infinity). */
		double max_timestep;
		/** Relative tolerance of the solver. */
		double rel_tol;
		/** Absolute tolerance of the solver. */
		double abs_tol;
	};

	/** A structure for the result of an initial value problem */
	struct ivp_result {
		typedef enum {
			COMPLETED, INTERRUPTED, ERROR
		} status_type;
		status_type status;
		scalar_type stop_time;
		state stop_state;
		trajectory<scalar_type> traj;

		template<typename result_type> typename ode_solver<result_type>::ivp_result convert_to() const {
			typename ode_solver<result_type>::ivp_result res;
			res.status
					= typename ode_solver<result_type>::ivp_result::status_type(
							status);
			res.stop_time = convert_element<result_type> (stop_time);
			res.stop_state = stop_state.template convert_to<result_type> ();
			res.traj = traj.template convert_to<result_type> ();
			return res;
		}
		;
	};

	/** A structure for the result of a root finding problem */
	struct rootf_result: ivp_result {
		rootf_result() {
		}
		;
		rootf_result(const ivp_result& res) :
			ivp_result(res) {
		}
		;
		typedef enum {
			FROM_BELOW, FROM_ABOVE, NO_ROOT
		} root_status;
		/** A row for every root found. */
		typedef math::matrix<root_status> root_matrix;
		root_matrix root_info;

		template<typename result_type> typename ode_solver<result_type>::rootf_result convert_to() const {
			typedef typename ode_solver<result_type>::rootf_result res_type;
			res_type res(this->ivp_result::template convert_to<result_type>());
			// This is a silly type conversion because gcc doesn't figure out it's actually
			// the same matrix type. Oh well...
			res.root_info = typename res_type::root_matrix(
					root_info.template convert_to<
							typename res_type::root_status> ());
			return res;
		}
		;
	};

	/** Initialize the initial value problem.
	 *
	 * To be run before solve_ivp or solve_roots is called.
	 *
	 * The ivp solves the ode given by f starting from state x0
	 * at time t0 up to time tfinal.
	 *
	 * dyn is a state_functor, which provides a function f(x):X->X
	 * for computing the derivative in state x.
	 *
	 * Calls init_ivp_impl to initialize derived classes.
	 * */
	virtual void init_ivp(const state& x0, scalar_type t0,
			const state_functor& f, ode_parameters par = ode_parameters());

	/** Solve the currently defined initial value problem up to time tfinal.
	 * May return just the final state, even when tfinal>max_timestep.
	 */
	virtual ivp_result step_ivp(scalar_type tfinal) = 0;

	/** Solve the currently defined initial value problem up to time tfinal.
	 *
	 * Returns a trajectory with all states computed up to tfinal,
	 * including the initial states.
	 */
	virtual ivp_result solve_ivp(scalar_type tfinal);

	/** Initialize the root finding problem.
	 *
	 * To be run before solve_ivp or solve_roots is called.
	 *
	 * The ivp solves the ode given by f starting from state x0
	 * at time t0 up to time tfinal, and during the solving
	 * checks for roots of a function g(x):X->Y on an arbitrary domain Y.
	 *
	 * dyn is a state_functor, which provides a function f(x):X->X
	 * for computing the derivative in state x.
	 *
	 * Calls init_rootf_impl to initialize derived classes.
	 * (init_ivp_impl is not called).
	 * */
	virtual void init_rootf(const state& x0, scalar_type t0,
			const state_functor& f, const state_vector_functor& g, ode_parameters par =
					ode_parameters());

	/**
	 *
	 * Specify a functor for Jacobian matrix computation
	 */

	virtual void specify_jacobian(const state_matrix_functor & j) = 0;


	/** Solve the currently defined root finding problem up to time tfinal,
	 * stopping as soon as a root is found.
	 *
	 * May return just the final state, even when tfinal>max_timestep.
	 * Returns at most one root.
	 */
	virtual rootf_result step_rootf(scalar_type tfinal) = 0;

	/** Solve the currently defined root finding problem up to time tfinal.
	 * Returns all roots found from t0 (not including) up to tfinal (including).
	 *
	 * The ith state in the returned trajectory corresponds to the
	 * ith root in root_info.
	 *
	 * @note To chain several calls to solve_rootf, use stop_state.
	 */
	virtual rootf_result solve_rootf(scalar_type tfinal);

	/** Solve the currently defined root finding problem up to time tfinal, or stops
	 *  at soon as a root is found.
	 *
	 * @param 	tfinal, stop time if no root is found prior.
	 * @return 	The trajectory with all computed states, not including initial state.
	 * 			if a root is found, it is the last one (stop_state) and the root's
	 * 			information are added to root_info.
	 *
	 */
	virtual rootf_result solve_firstrootf(scalar_type tfinal);

	/** Get the current time of integration of the solver*/
	virtual scalar_type get_solver_current_time() const = 0;

	/** Get a string name for the solver. */
	virtual std::string get_solver_name() const = 0;

protected:
	/** Further initialization by the solver implementation. */
	virtual void init_ivp_impl() = 0;
	/** Further initialization by the solver implementation. */
	virtual void init_rootf_impl() = 0;

	state my_x0;
	scalar_type my_t0;
	const state_functor* my_f;
	const state_vector_functor* my_g;
	const state_matrix_functor * my_jac;
	ode_parameters my_odepar;
};

}
}

/** Stream output for lp_result. */
template<typename scalar_type> std::ostream& print(std::ostream& os,
		const typename math::ode::ode_solver<scalar_type>::ivp_result& res) {
	if (res.status == math::ode::ode_solver<scalar_type>::ivp_result::ERROR) {
		os << "ERROR" << std::endl;
	}
	os << res.traj;
	return os;
}
;

#include "ode_solver.hpp"

#endif /* ODE_SOLVER_H_ */
