/*
 * traj_simu.h
 *
 *  Created on: 10 mars 2011
 *      Author: christophe
 */

#ifndef TRAJ_SIMU_H_
#define TRAJ_SIMU_H_

#include "math/vdom/vdom_vector.h"
#include "math/vdom/affine_map.h"
#include "math/vdom/trajectory.h"
#include "math/vdom/lin_constraint_system.h"
#include "math/ode_solving/ode_solver_cvode.h"
#include "math/numeric/interval.h"
#include "math/numeric/approx_comparator.h"
#include "math/numeric/comp.h"
#include "core/hybrid_automata/location_id.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "math/vdom/state_vector_functor_utility.h"

template<typename scalar_type>
class traj_simu {
public:

	// typedefs for the algorithm
	typedef math::numeric::approx_comparator<scalar_type> approx_comparator;

	typedef math::vdom_vector<scalar_type> state;

	typedef math::numeric::interval<scalar_type> interval;

	typedef math::trajectory<scalar_type> trajectory;

	typedef math::affine_map<scalar_type> affine_map;

	typedef math::lin_constraint_system<scalar_type> lin_constraint_system;

	typedef  math::lin_constraint_state_functor<scalar_type>
														lin_constraint_state_functor;

	typedef math::affine_map_vector_state_vector_functor<scalar_type>
														affine_map_vector_state_vector_functor;

	typedef typename math::state_functor<scalar_type> state_functor;
	typedef typename math::state_matrix_functor<scalar_type> state_matrix_functor;

	typedef typename lin_constraint_state_functor::lin_constraints lin_constraints;
	typedef typename lin_constraint_state_functor::lin_constraints_vector lin_constraints_vector;

	typedef math::ode::ode_solver_cvode<scalar_type> ode_solver;

	/**
	 * A jump interval is a time interval, with indexes corresponding
	 *  to positions in the corresponding trajectory
	 */
	typedef struct jump_interval{
		interval interv;
		unsigned int index_low;
		unsigned int index_high;
		//friend std::ostream& operator<< (std::ostream& os, const traj_simu<scalar_type>::jump_interval & intv);
	} jump_interval;

	typedef std::list<jump_interval> interval_list;

	/**
	 *  Structure used as the return of state_jump_intervals.
	 */
	typedef struct state_jump_results{
		bool horizont_reached;
		std::vector<interval_list> ji;
		trajectory traj;
	} state_jump_results;


	/**
	 * Basically a vector of collection of (timed) states - used to represent the future roots of the next iteration of the algorithm.
	 */
	typedef std::vector<trajectory> simulation_results;

	/**
	 * Structure used to store the result of a whole simulation :  trajectory + the new set of roots selected for each guard ( simulation
	 * will begin again in the locations of each transition's target with the selected roots for this transition)
	 */
	typedef struct roots_and_trajectory{
		trajectory traj;
		simulation_results roots;

	} roots_and_trajectory;

	/**
	 * Parameters used by the simulation algorithm, currently only global and
	 * local time horizon.
	 */
	typedef struct traj_simu_parameters{
		scalar_type global_horizont;
		scalar_type local_horizont;

		traj_simu_parameters(){
			global_horizont = scalar_type(-1);
			local_horizont = scalar_type(-1);
		}
		traj_simu_parameters(scalar_type gh, scalar_type lh){
			global_horizont = gh;
			local_horizont = lh;
		}
	} traj_simu_parameters;

	/**
	 *
	 * Constructor that applies default parameters, e.g. no time horizons and a time step of 0.1
	 * Error bounds remains as previously set for approx_comparator and ode_solver.
	 */
	traj_simu() :
		par(), params()
	{

		par.max_timestep=0.1;
		/*
		approx_comparator::set_rel_error_bound(1e-8);
		approx_comparator::set_abs_error_bound(1e-10);
		par.rel_tol = approx_comparator::get_rel_error_bound()*0.0001;
		par.abs_tol = approx_comparator::get_abs_error_bound()*0.0001;
		*/

//		LOGGER_OS(DEBUG7,"postd_simulation") << "DEBUG traj_simu : ODE : rel_tol = "<< par.rel_tol <<" abs_tol = " << par.abs_tol;
//		LOGGER_OS(DEBUG7,"postd_simulation") << "DEBUG traj_simu : program : rel_tol = "<< approx_comparator::get_rel_error_bound() <<" abs_tol = " << approx_comparator::get_abs_error_bound();
	};


	/**
	 *
	 * Constructor taking user-defined parameters for the solver and the simulation
	 * (such as time step, tolerance values, time horizons etc...)
	 * For the solver parameters, see corresponding documentation
	 */
	traj_simu(typename ode_solver::ode_parameters par, traj_simu_parameters params ){

		this->par = par;
		this->params = params;
		/*
		if(this->par.max_timestep>0.1 || this->par.max_timestep<0) {
				this->par.max_timestep = 0.01;
				printf("traj_simu::constructor reset granularity to 0.01\n");
		}

		approx_comparator::set_rel_error_bound(1e-8);
		approx_comparator::set_abs_error_bound(1e-10);
		this->par.rel_tol = approx_comparator::get_rel_error_bound()*0.0001;
		this->par.abs_tol = approx_comparator::get_abs_error_bound()*0.0001;
		 */

//		LOGGER_OS(DEBUG7,"postd_simulation") << "DEBUG traj_simu : ODE : rel_tol = "<< par.rel_tol <<" abs_tol = " << par.abs_tol;
//		LOGGER_OS(DEBUG7,"postd_simulation") << "DEBUG traj_simu : program : rel_tol = "<< approx_comparator::get_rel_error_bound() <<" abs_tol = " << approx_comparator::get_abs_error_bound();
	};

	 ~traj_simu() {};

	 /**
	 	 * A function that takes an initial time and state, and with the help of the derivative x'=f(x,t) integrates the trajectory up to tf, or stops as soon
	 	 * as the invariant becomes unsatisfied. Additionally, it checks for roots of the guards, and returns the intervals in which the different guards are
	 	 * satisfied. An empty interval list on the guard i means that the guard i is never satisfied during the trajectory from tinit to tfinal (that is <= tf
	 	 * obviously).
	 	 *
	 	 * This variant consults the solver for informations on the directions on roots.
	 	 * Both g(x) and g'(x) are checked for roots, in order to grant better detection of atomic points (though this detection is not guaranteed).
	 	 *
	 	 * @param f_i functor calculating x_dot = f(x)
	 	 * @param Jac_i functor calculating the Jacobian
	 	 * @param xinit the intial value for x
	 	 * @param tinit the initial time e.g. at init xinit = x(tinit)
	 	 * @param tmax  the upper bound of the integration's time
	 	 * @param guards a vector containing all guards in linear constraint form, of size n
	 	 * @param invariant the invariant in linear constriant form
	 	 * @return 	a structure containing the integrated trajectory, and a vector ji, where ji(i) is a vector containing the intervals during which the i-th guard
	 	 * 			is satisfied
	 	 */
	 	state_jump_results compute_urgent_traj(const state_functor & f_i,  const state_matrix_functor & Jac_i,  const state & xinit,
	 			scalar_type tinit, scalar_type tmax, std::vector<typename lin_constraint_system::const_ptr> & guards,
	 			typename lin_constraint_system::const_ptr invariant ) const;


	 /**
	 	 * A function that takes an initial time and state, and with the help of the derivative x'=f(x,t) integrates the trajectory up to tf, or stops as soon
	 	 * as the invariant becomes unsatisfied. Additionally, it checks for roots of the guards, and returns the intervals in which the different guards are
	 	 * satisfied. An empty interval list on the guard i means that the guard i is never satisfied during the trajectory from tinit to tfinal (that is <= tf
	 	 * obviously).
	 	 *
	 	 * This variant consults the solver for informations on the directions on roots.
	 	 * Both g(x) and g'(x) are checked for roots, in order to grant better detection of atomic points (though this detection is not guaranteed).
	 	 *
	 	 * @param f_i functor calculating x_dot = f(x)
	 	 * @param Jac_i functor calculating the Jacobian
	 	 * @param xinit the intial value for x
	 	 * @param tinit the initial time e.g. at init xinit = x(tinit)
	 	 * @param tmax  the upper bound of the integration's time
	 	 * @param guards a vector containing all guards in linear constraint form, of size n
	 	 * @param invariant the invariant in linear constriant form
	 	 * @return 	a structure containing the integrated trajectory, and a vector ji, where ji(i) is a vector containing the intervals during which the i-th guard
	 	 * 			is satisfied
	 	 */
	 	state_jump_results state_jump_intervals_solver(const state_functor & f_i,  const state_matrix_functor & Jac_i,  const state & xinit,
	 			scalar_type tinit, scalar_type tmax, std::vector<typename lin_constraint_system::const_ptr> & guards,
	 			typename lin_constraint_system::const_ptr invariant ) const;

	/**
	 * See  state_jump_intervals_solver, the major difference between the two functions is that the _solver variant
	 * consults the solver for informations on the directions on roots, whereas this one calculates them by evaluating
	 * the sign of g'(x).
	 * Detection of atomic roots is not guaranteed if g'(x) is zero-valued at the root.
	 *
	 *
	 * @param f_i functor calculating x_dot = f(x)
	 * @param Jac_i functor calculating the Jacobian
	 * @param xinit the intial value for x
	 * @param tinit the initial time e.g. at init xinit = x(tinit)
	 * @param tmax  the upper bound of the integration's time
	 * @param guards a vector containing all guards in linear constraint form, of size n
	 * @param invariant the invariant in linear constraint form
	 * @return 	a structure containing the integrated trajectory, and a vector ji, where ji(i) is a vector containing the intervals during which the i-th guard
	 * 			is satisfied
	 */
	state_jump_results state_jump_intervals_derivative(const state_functor & f_i,  const state_matrix_functor & Jac_i,  const state & xinit,
			scalar_type tinit, scalar_type tmax, std::vector<typename lin_constraint_system::const_ptr> & guards,
			typename lin_constraint_system::const_ptr invariant ) const;

	/**
	 * A function that selects next roots by discretizing the jump intervals found by a former call to one of the state_jump_intervals_...
	 * This is a debug/simple version of the function : it only
	 * select 3 points : lower + upper bound (t) and a  sampled point in the middle
	 *
	 * @param res : results returned by a former call to state_jump_interval_xxx
	 * @return  r the list of new possible states with possible values of x
	 * 			r.traj  = 	trajectory computed by the ODE solver
	 * 			r.roots =	selected roots for the next iteration. WARNING : This roots have not yet been imaged thorugh the transition, see
	 * 						image_next_traj_states
	 */
	roots_and_trajectory select_next_traj_states( const state_jump_results & res) const ;


	/**
	 *
	 * @param res : simulation_results obtained via a former call to select_next_traj_states
	 * @param tr_assign : a vector of functors that is of the same size as res.
	 * @return a vector of trajectories of the same size as res (if res[i] is empty, the i-th trajectory of the result is empty
	 */
	simulation_results image_next_traj_states( const simulation_results & res, const  std::vector<typename state_functor::const_ptr> & tr_assign) const ;

	/**
	 *	Calculates the image of a state through a function (corresponding here to the assignment of a transition)
	 * @param s : a state on a domain d
	 * @param f_assign : a functor capable of calculating s' == f_assign(s) ( dom(s') should be d )
	 * @return s' == f_assign(s)
	 */
	state state_image_through_transition(const state & s, const state_functor & f_assign) const;


private:
	typename ode_solver::ode_parameters par;
	traj_simu_parameters params;
};

/**
 * A state to state functor that reorder the domain of the result according to the input state's domain
 */
template<typename scalar_type>
class reorder_state_functor : public math::state_functor<scalar_type> {
public:
	typedef typename math::state_functor<scalar_type>::state state;

	reorder_state_functor(const math::state_functor<scalar_type> & f) : my_reorderf(f){
	}

	virtual ~reorder_state_functor(){
	}

	virtual state map(const state& x) const{
		state temp =my_reorderf.map(x);
		temp.reorder(x.domain());
		return temp;
	}

private:
	 const math::state_functor<scalar_type> & my_reorderf;
};

#include "traj_simu.hpp"



/**
 *
 * This code is not working, because jump interval is a class encapsulated in a
 * template class.
 * This code should be instantiated according to the scalar_type used
 * See function underneath
 */
template<typename scalar_type>
std::ostream& operator<<(std::ostream& os, const  typename traj_simu<
		scalar_type>::jump_interval & intv) {
	os << "[" << intv.interv.lower() << "," << intv.interv.upper() << "]";
   return os;
}

/**
 * Ostream operator for jump_interval, instantiated for double
 */
inline std::ostream& operator<<(std::ostream& os, const traj_simu<
		double>::jump_interval & intv) {
	os << "[" << intv.interv.lower() << "," << intv.interv.upper() << "]";
   return os;
}




#endif /* TRAJ_SIMU_H_ */
