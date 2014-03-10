/*
 * ode_solver.hpp
 *
 *  Created on: Oct 11, 2010
 *      Author: frehse
 */

#ifndef ODE_SOLVER_HPP_
#define ODE_SOLVER_HPP_

#include "ode_solver.h"

namespace math {
namespace ode {

template<typename scalar_type>
ode_solver<scalar_type>::~ode_solver() {
}

template<typename scalar_type>
void ode_solver<scalar_type>::init_ivp(const state& x0, scalar_type t0,
		const state_functor& f, ode_parameters par) {
	my_x0 = x0;
	my_t0 = t0;
	my_f = &f;
	my_g = 0;
	my_odepar = par;
	init_ivp_impl();
}

template<typename scalar_type>
typename ode_solver<scalar_type>::ivp_result ode_solver<scalar_type>::solve_ivp(
		scalar_type tfinal) {

	ivp_result res, stepres;
	res.status = ivp_result::COMPLETED;

	res.traj = trajectory<scalar_type> (my_t0, my_x0);
	scalar_type t =this->get_solver_current_time();
	scalar_type tf;
	while (t < tfinal && res.status == ivp_result::COMPLETED) {
		if (my_odepar.max_timestep > 0.0)
			tf = t + scalar_type(my_odepar.max_timestep);
		else
			tf = tfinal;
		/* Solve for one timestep */
		stepres = step_ivp(tf);
		/* Join the new state with result so far. */
		res.status = stepres.status;
		res.stop_time = stepres.stop_time;
		res.traj.insert(stepres.traj);
		t = stepres.stop_time;
	}
	res.stop_state = stepres.stop_state;
	return res;
}

template<typename scalar_type>
void ode_solver<scalar_type>::init_rootf(const state& x0, scalar_type t0,
		const state_functor& f, const state_vector_functor& g, ode_parameters par) {
	my_x0 = x0;
	my_t0 = t0;
	my_f = &f;
	my_g = &g;
	my_odepar = par;
	init_rootf_impl();
}

template<typename scalar_type>
typename ode_solver<scalar_type>::rootf_result ode_solver<scalar_type>::solve_rootf(
		scalar_type tfinal) {

	rootf_result res, stepres;
	res.status = rootf_result::COMPLETED;

	res.traj = trajectory<scalar_type> (my_x0.domain());
	scalar_type t = this->get_solver_current_time();
	scalar_type tf;
	bool last_state = false;
	while (!last_state && res.status != ivp_result::ERROR) {
		if (my_odepar.max_timestep > 0.0)
			tf = t + scalar_type(my_odepar.max_timestep);
		else
			tf = tfinal;
		/* Solve for one timestep */
		stepres = step_rootf(tf);
		/* Join the new state with result so far. */
		res.status = stepres.status;
		res.stop_time = stepres.stop_time;
		t = stepres.stop_time;
//std::cout << stepres.stop_time << " with " << stepres.status << std::endl;

		/** It's the last state if stop_time is sufficiently close to tfinal */
		last_state = (t >= tfinal*scalar_type(1.0-my_odepar.rel_tol)-scalar_type(my_odepar.abs_tol));

		/** Include stop state if it's a root. */
		if (stepres.status == rootf_result::INTERRUPTED) {
			unsigned int last_row = stepres.traj.size() - 1;
			res.traj.insert(stepres.traj, last_row);
			if (stepres.root_info.size1() > 1)
				throw std::runtime_error(
						"step_rootf returned more than one root");
			res.root_info.append_row(stepres.root_info.vector_from_row(0));
		}
	}
	res.stop_state = stepres.stop_state;
	return res;
}


template<typename scalar_type>
typename ode_solver<scalar_type>::rootf_result ode_solver<scalar_type>::solve_firstrootf(
		scalar_type tfinal) {

	rootf_result res, stepres;
	res.status = rootf_result::COMPLETED;

	res.traj = trajectory<scalar_type> (my_x0.domain());
	scalar_type t = this->get_solver_current_time();
	scalar_type tf;
	int sign = ((tfinal - t)>0)?1:-1;
	bool last_state = false;
	unsigned int pass = 0;
	while (!last_state && stepres.status != ivp_result::ERROR) {

		pass++;

		if (my_odepar.max_timestep > 0.0){
			tf = t + scalar_type(sign)*scalar_type(my_odepar.max_timestep);
			if( sign*tf>=sign*tfinal) tf = tfinal;
		}
		else
			tf = tfinal;



		/* Solve for one timestep */
		stepres = step_rootf(tf);
		/* Join the new state with result so far. */
		t = stepres.stop_time;
		res.traj.insert(stepres.traj);
		//std::cout << stepres.stop_time << " with " << stepres.status << std::endl;




		/** It's the last state if stop_time is sufficiently close to tfinal */
		if(tfinal>=scalar_type(0) )
		{
			if(sign>0)
				last_state = ( t >= (tfinal*scalar_type(1.0-my_odepar.rel_tol)-scalar_type(my_odepar.abs_tol))) ;
			else
				last_state = ( t<= (tfinal*scalar_type(1.0+my_odepar.rel_tol)+scalar_type(my_odepar.abs_tol))) ;

		}
		else
		{
			if(sign>0)
				last_state = ( t >= (tfinal*scalar_type(1.0+my_odepar.rel_tol)-scalar_type(my_odepar.abs_tol))) ;
			else
				last_state = ( t<= (tfinal*scalar_type(1.0-my_odepar.rel_tol)+scalar_type(my_odepar.abs_tol))) ;
		}
		/** if it's a root, stop ode solving and add root infos. */
		if (stepres.status == rootf_result::INTERRUPTED) {

			res.root_info.append_row(stepres.root_info.vector_from_row(0));
			last_state=true;

			if (stepres.root_info.size1() > 1)
				throw std::runtime_error(
						"step_rootf returned more than one root");

			}

	}
	res.status = stepres.status;
	res.stop_time = stepres.stop_time;
	res.stop_state = stepres.stop_state;
	//my_x0=res.stop_state;
	//my_t0=res.stop_time;
	//printf("solve_firstrootf : stop_time = %e \n",res.stop_time);
	return res;
}

/*****************************************************************/

}
}

#endif /* ODE_SOLVER_HPP_ */
