/*
 * traj_simu.cpp
 *
 *  Created on: 10 mars 2011
 *      Author: christophe
 */

#ifndef TRAJ_SIMU_HPP_
#define TRAJ_SIMU_HPP_

#include "traj_simu.h"

#include "math/vdom/state_functor_utility.h"
#include "math/vdom/lin_constraint_system_utility.h"
#include "math/vdom/vdom_matrix.h"
#include "math/tribool.h"

//debug
#include "utility/stl_helper_functions.h"
#include "utility/basic_warning.h"
#include <iostream>
#include <sstream>

using namespace math;
using namespace numeric;


//TODO finish this
template<typename scalar_type>
typename traj_simu<scalar_type>::state_jump_results traj_simu<scalar_type>::compute_urgent_traj(
		const state_functor & f, const state_matrix_functor & Jac,
		const state & xinit, scalar_type tinit, scalar_type tmax,
		std::vector<typename lin_constraint_system::const_ptr> & guards,
		typename lin_constraint_system::const_ptr invariant) const {

	//used to store intermediate IVP and Root results returned by the solver
	typename ode_solver::rootf_result interm_result;

  // typedef for root finding
	typedef typename math::ode::ode_solver<scalar_type>::rootf_result::root_matrix root_matrix;
	typedef typename ode_solver::rootf_result::root_status root_status;

	//solver constants for rootfinding

	root_status rsNO_ROOT = ode_solver::rootf_result::NO_ROOT;
	root_status rsFROM_BELOW = ode_solver::rootf_result::FROM_BELOW;
	root_status rsFROM_ABOVE = ode_solver::rootf_result::FROM_ABOVE;

	unsigned int gn = guards.size(); // number of guards

	// variables used to store temporary/intermediate results results
	jump_interval temp;
	int temp_satisfied;
	scalar_type xp; //scalar used for initialisation
	state xinitprim;
	scalar_type zero(0);

	//the result of this function
	state_jump_results res;

	// variables used to determine which guard is satisfied or not, and if it is, since when, or if it's not, how many
	// constraints are not yet satisfied
	scalar_type tdeb_k[gn]; //if jump_ok[i] then tdeb_k[i] time the root occurred
	unsigned int tdeb_index[gn];
	int remaining_to_satisf[gn];
	bool jump_ok[gn]; //is the jump OK at current time ?

	//variables used to pass constraints of guards (and invariant) to the solver and
	//memorize indexes in order to assign results produced by the solver
	//to the corresponding guard ( or the invariant)
	lin_constraints_vector z(gn + 1); //current sought roots
	int indices[2 * (gn + 1)];

	//loop variable (iterators)
	typename lin_constraint_system::const_iterator it;
	typename lin_constraint_system::const_iterator git;
	typename lin_constraint_system::const_iterator oit;

	//boolean for guards
	bool guard_satisfied = true;
	bool guard_atomic = false;

	// an instance of the SOLVER --> needs an initialisation (performed later)
	ode_solver solver;

	///checking invariant + which guard is sastified and which is'nt.
	///guards are assumed of the form a*x+b<=0

	/// Filling z and indices
	/// with invariant and guards constraints

	// to bring lin_constraint_system to domain of solver
	typename lin_constraint_system::ptr reordered_invariant =
			typename lin_constraint_system::ptr(
					new lin_constraint_system(*invariant));
	// @todo choose a domain for the solver and stick with it (use it everywhere)
	reordered_invariant->reorder(xinit.domain());
	invariant = reordered_invariant;

	for (int i = 0; i < gn; i++) {
		typename lin_constraint_system::ptr reordered_guard =
				typename lin_constraint_system::ptr(
						new lin_constraint_system(*guards[i]));
		// @todo choose a domain for the solver and stick with it (use it everywhere)
		reordered_guard->reorder(xinit.domain());
		guards[i] = reordered_guard;
	}

	z[0].begin = invariant->begin();
	z[0].end = invariant->end();
	indices[0] = 0;
	indices[1] = 2 * invariant->size();

	for (int i = 0; i < gn; i++) {
		LOGGER_OS(DEBUG7,"compute_urgent_traj")
			<< "guard " << i << ":\n" << guards[i];
		//setting z for this guard
		z[i + 1].begin = guards[i]->begin();
		z[i + 1].end = guards[i]->end();

		indices[(i + 1) * 2] = indices[(i + 1) * 2 - 1];
		indices[(i + 1) * 2 + 1] = indices[(i + 1) * 2] + 2 * guards[i]->size();
	}

	lin_constraint_state_functor z_functor(z,
			(const math::state_functor<scalar_type> *) &f, indices[2 * gn + 1]);

	/// solver initialisation
	solver.init_rootf(xinit, tinit, f, z_functor, par);
	solver.specify_jacobian(Jac);

	///do solver adjustment :
	/// a little backward then forward integration up to t = tinit in
	/// order to gather information on roots of invariants and guards at t = tinit
	/// (to detect and handle properly roots at the beginning of the algorithm)
	interm_result = solver.solve_firstrootf(tinit - 1);
	while (!approx_comparator::is_definitely_strictly_smaller(
			interm_result.stop_time, tinit)) {
		interm_result = solver.solve_firstrootf(tinit - 1);
	}
	solver.init_rootf(interm_result.stop_state, interm_result.stop_time, f,
			z_functor, par);
	while (approx_comparator::is_definitely_strictly_smaller(
			interm_result.stop_time, tinit)) {
		interm_result = solver.solve_firstrootf(tinit + par.max_timestep);
	}

	/// if solver has not stopped at t = tinit, there is no root at the beginning of the integration.
	if (approx_comparator::is_definitely_strictly_larger(
			interm_result.stop_time, tinit) && interm_result.root_info.size1()> 0) {
		interm_result.root_info = root_matrix(1, indices[gn * 2 + 1], rsNO_ROOT);
	}

	///for verification of invariant
	/// sat is true if the invariant is true, atom is true is the invariant is atomic, e.g. only satisfied at t=tinit
	bool invariant_satisfied = true;
	bool invariant_atomic = false;

	/// root_info matrix should not be empty in the next loops
	if (interm_result.root_info.size1() == 0)
		interm_result.root_info
				= root_matrix(1, indices[gn * 2 + 1], rsNO_ROOT);

	/// evaluating invariant and then guards values and directions (if roots)  at t=tinit
	/// Goal is to detect if a guard/ the invariant is satisfied and if it is the case, verify
	/// if this situation is atomic, e.g. one of the constraint is rooted and goes up, which
	/// means that the constraint set will not be satisfied after t = current time
	/// The same sort of check is performed in the main iteration each time the solver stops with a root.
	it = invariant->begin();
	for (int i = 0; (2 * i) < indices[1]; i++) {
		LOGGER_OS(DEBUG7,"compute_urgent_traj")
			<< "initial check on invariant " << i << ":" << *it << " : ";
		// Check whether x violates the constraint or is on the border
		if (interm_result.root_info(0, 2 * i) == rsFROM_BELOW) {
			invariant_atomic = true;
		} else if (interm_result.root_info(0, 2 * i) == rsNO_ROOT) {
			xp = lin_constraint_evaluation(*it, xinit);
			if (approx_comparator::is_definitely_strictly_pos(xp)) {
				invariant_satisfied = false;
			} else if (!approx_comparator::is_definitely_strictly_neg(xp)) {
				// The invariant constraint is on the border.
				// Check the derivative to see the whether the invariant function is increasing or decreasing.
				if (interm_result.root_info(0, 2 * i + 1) == rsFROM_BELOW) {
					invariant_atomic = true;
				}
			}
		}
		++it;
	}

	/// If the invariant is not satisfied at the beginning, no trajectory should be computed and no guards taken
	/// (state is invalid !!!)
	if (!invariant_satisfied) {
		LOGGER_OS(DEBUG5,"compute_urgent_traj")
					<< "Invariant is not satisfied at the beginning of the simulation !!!";
		res.traj = trajectory(xinit.domain());
		res.ji = std::vector<interval_list>(gn);
		res.horizont_reached = false;
		return res;
	}

	/// If the invariant is not atomic at the beginning, no further trajectory will be computed but guards may be taken.
	/// See condition of main loop
	if (invariant_atomic)
		LOGGER_OS(DEBUG5,"compute_urgent_traj")
			<< "Inv atomic at beginning !!!";

	/// now evaluating if guards are satisfied at beginning.
	/// the atomic et satisfied booleans play the same role as atom & satisf play for the infariant
	res.traj = trajectory(xinit.domain(), tinit, xinit.get_vector());
	res.ji = std::vector<interval_list>(gn);

	for (int i = 0; i < gn; i++) {
		it = guards[i]->begin();
		remaining_to_satisf[i] = 0;
		guard_atomic = false;
		guard_satisfied = true;
		for (int j = 0; j < guards[i]->size(); j++) {
			LOGGER_OS(DEBUG7,"compute_urgent_traj")
				<< "initial check on guard " << i << ":" << *it << " : ";

			if (interm_result.root_info(0, indices[2 * (i + 1)] + 2 * j)
					== rsFROM_BELOW) {
				LOGGER_OS(DEBUG7,"compute_urgent_traj")
					<< "atomic" << "\n";
				guard_atomic = true;
				remaining_to_satisf[i]++;

			} else if (interm_result.root_info(0, indices[2 * (i + 1)] + 2 * j)
					== rsNO_ROOT) {
				xp = lin_constraint_evaluation(*it, xinit);
				scalar_type inh_coeff = it->get_canonic_inh_coeff();
				if (definitely(is_GT(xp - inh_coeff, -inh_coeff))) {
					LOGGER_OS(DEBUG7,"compute_urgent_traj")
						<< "violated by " << xp << "\n";

					guard_satisfied = false;
					remaining_to_satisf[i]++;
				} else if (!definitely(is_LT(xp - inh_coeff, -inh_coeff))) {
					LOGGER_OS(DEBUG7,"compute_urgent_traj")
						<< "zero at " << xp << "\n";

					// Check the derivative to see the whether the guard function is increasing or decreasing.
					// The derivative is f(x). If a^Tf(x)>0, then the function is increasing.
					if (interm_result.root_info(0,
							indices[2 * (i + 1)] + 2 * j + 1) == rsFROM_BELOW) {
						LOGGER_OS(DEBUG7,"compute_urgent_traj")
							<< "from below, so atomic " << "\n";
						guard_atomic = true;
						// The constraint is violated for all t'>t, so
						// it needs to be added to the constraints remaining to satisfy
						remaining_to_satisf[i]++;

					}
				}
			}

			++it;
		}
		if (guard_satisfied && guard_atomic) {
			jump_ok[i] = false;
			temp.interv.set_lower(tinit);
			temp.interv.set_upper(tinit);
			temp.index_low = 0;
			temp.index_high = 0;
			res.ji[i].push_back(temp);
		} else if (guard_satisfied) {
			jump_ok[i] = true;
			tdeb_k[i] = tinit;
			tdeb_index[i] = 0;
		} else
			jump_ok[i] = false;
	}

	solver.init_rootf(xinit, tinit, f, z_functor, par); // init solver
	solver.specify_jacobian(Jac);

	//algorithm main loop
	res.horizont_reached = false;
	bool inv = !invariant_atomic;
	scalar_type t = tinit;

	scalar_type tmax_next;
	scalar_type tprec;
	// inv && tmax>=t
	while (inv && (tmax < 0 || !maybe(approx_comparator::is_LE(tmax, t)))) {

		//TODO find a better solution if tmax < 0
		tmax_next = (tmax < 0) ? (2 * t + 1) : tmax;

		interm_result = solver.solve_firstrootf(tmax_next);
		res.traj.insert(interm_result.traj);
		tprec = t;
		t = interm_result.stop_time;

		if (!(interm_result.root_info.size1() == 0)
				&& approx_comparator::is_definitely_strictly_larger(t, tprec)) {

			LOGGER_OS(DEBUG7,"compute_urgent_traj")
				<< "Root infos : NO_ROOT=" << (int) rsNO_ROOT
						<< ", FROM_ABOVE=" << (int) rsFROM_ABOVE
						<< ", FROM_BELOW=" << (int) rsFROM_BELOW << std::endl
						<< "root matrix returned by solver "
						<< interm_result.root_info << "and g values : "
						<< z_functor.map(interm_result.stop_state);

			it = invariant->begin();
			for (int i = 0; (2 * i) < indices[1]; i++) {
				xp = lin_constraint_evaluation(*it, interm_result.stop_state);
				if (interm_result.root_info(0, 2 * i) == rsFROM_BELOW) {
					inv = false;
					LOGGER_OS(DEBUG7,"compute_urgent_traj")
						<< "invariant rooted, index = " << 2 * i
								<< "with value " << xp;
				} else LOGGER_OS(DEBUG7,"compute_urgent_traj")
					<< "invariant not rooted, index = " << 2 * i
							<< "with value " << xp;
				++it;
			}

			for (int i = 0; i < gn; i++) //for all guards
			{

				it = guards[i]->begin();
				guard_atomic = false;
				temp_satisfied = 0;
				for (int j = 0; j < guards[i]->size(); j++) {
					LOGGER_OS(DEBUG7,"compute_urgent_traj")
						<< "check on guard " << i << ", index=" << (indices[2
								* (i + 1)] + 2 * j) << " :" << *it << " : ";

					// Note: The solver seems to have quite unreliable root info
					// (timesteps too big?)
					// So let's try using the derivative

					bool below = false;
					bool above = false;

					below = interm_result.root_info(0,
							indices[2 * (i + 1)] + 2 * j) == rsFROM_BELOW;
					above = interm_result.root_info(0,
							indices[2 * (i + 1)] + 2 * j) == rsFROM_ABOVE;

					if (below) {
						LOGGER_OS(DEBUG7,"compute_urgent_traj")
							<< "from below, so atomic";

						guard_atomic = true;
						remaining_to_satisf[i]++;
						temp_satisfied++;

					} else if (above) {
						remaining_to_satisf[i]--;

						LOGGER_OS(DEBUG7,"compute_urgent_traj")
							<< "from above, remaining:"
									<< remaining_to_satisf[i];
					} else {

						xp = lin_constraint_evaluation(*it,
								interm_result.stop_state);
						scalar_type inh_coeff = it->get_canonic_inh_coeff();
						double dxp = lin_constraint_evaluation(*it,
								f.map(interm_result.stop_state));
						LOGGER_OS(DEBUG7,"compute_urgent_traj")
							<< "solver says not root, guard value " << xp
									<< ", checking derivative dxp=" << dxp;

						if (interm_result.root_info(0,
								indices[2 * (i + 1)] + 2 * j + 1)
								== rsFROM_BELOW && is_MEQ(xp - inh_coeff,
								-inh_coeff)) {
							LOGGER_OS(DEBUG7,"compute_urgent_traj")
								<< "is satisfied";
							temp_satisfied++;
						}

					}

					++it;
				}

				if (jump_ok[i] && remaining_to_satisf[i] != 0) {

					temp.interv.set_lower(tdeb_k[i]);
					temp.interv.set_upper(t);
					temp.index_low = tdeb_index[i];
					temp.index_high = res.traj.size() - 1;
					res.ji[i].push_back(temp);

					jump_ok[i] = false;

				} else if (!jump_ok[i] && remaining_to_satisf[i] == 0) {
					jump_ok[i] = true;
					tdeb_k[i] = t;
					tdeb_index[i] = res.traj.size() - 1;

				} else if (!jump_ok[i] && remaining_to_satisf[i]
						- temp_satisfied == 0) {
					temp.interv.set_lower(t);
					temp.interv.set_upper(t);
					temp.index_low = res.traj.size() - 1;
					temp.index_high = temp.index_low;
					res.ji[i].push_back(temp);
				}

			}

		}
	}

	//if(t>=tmax) res.horizont_reached=true; -->tribool
	if ((!(tmax < 0)) && maybe(approx_comparator::is_LE(tmax, t)))
		res.horizont_reached = true;

	for (int i = 0; i < gn; i++) {
		if (jump_ok[i]) {
			temp.interv.set_lower(tdeb_k[i]);
			temp.interv.set_upper(t);
			temp.index_low = tdeb_index[i];
			temp.index_high = res.traj.size() - 1;
			res.ji[i].push_back(temp);
		}
	}

	return res;

}


//TODO finish this
template<typename scalar_type>
typename traj_simu<scalar_type>::state_jump_results traj_simu<scalar_type>::state_jump_intervals_solver(
		const state_functor & f, const state_matrix_functor & Jac,
		const state & xinit, scalar_type tinit, scalar_type tmax,
		std::vector<typename lin_constraint_system::const_ptr> & guards,
		typename lin_constraint_system::const_ptr invariant) const {

	//	reorder_state_functor<scalar_type> f = reorder_state_functor<scalar_type> (
	//			f_i);

	//used to store intermediate IVP and Root results returned by the solver
	typename ode_solver::rootf_result interm_result;

  // typedef for root finding
	typedef typename math::ode::ode_solver<scalar_type>::rootf_result::root_matrix root_matrix;
	typedef typename ode_solver::rootf_result::root_status root_status;

	//solver constants for rootfinding

	root_status rsNO_ROOT = ode_solver::rootf_result::NO_ROOT;
	root_status rsFROM_BELOW = ode_solver::rootf_result::FROM_BELOW;
	root_status rsFROM_ABOVE = ode_solver::rootf_result::FROM_ABOVE;

	unsigned int gn = guards.size(); // number of guards

	// variables used to store temporary/intermediate results results
	jump_interval temp;
	int temp_satisfied;
	scalar_type xp; //scalar used for initialisation
	state xinitprim;
	scalar_type zero(0);

	//the result of this function
	state_jump_results res;

	// variables used to determine which guard is satisfied or not, and if it is, since when, or if it's not, haw many
	// constraints are not yet satisfied
	scalar_type tdeb_k[gn]; //if jump_ok[i] then tdeb_k[i] time the root occurred
	unsigned int tdeb_index[gn];
	int remaining_to_satisf[gn];
	bool jump_ok[gn]; //is the jump OK at current time ?

	//variables used to pass constraints of guards (and invariant) to the solver and
	//memorize indexes in order to assign results produced by the solver
	//to the corresponding guard ( or the invariant)
	lin_constraints_vector z(gn + 1); //current sought roots
	int indices[2 * (gn + 1)];

	//loop variable (iterators)
	typename lin_constraint_system::const_iterator it;
	typename lin_constraint_system::const_iterator git;
	typename lin_constraint_system::const_iterator oit;

	//boolean for guards
	bool guard_satisfied = true;
	bool guard_atomic = false;

	// an instance of the SOLVER --> needs an initialisation (performed later)
	ode_solver solver;

	///checking invariant + which guard is sastified and which is'nt.
	///guards are assumed of the form a*x+b<=0

	/// Filling z and indices
	/// with invariant and guards constraints

	// to bring lin_constraint_system to domain of solver
	typename lin_constraint_system::ptr reordered_invariant =
			typename lin_constraint_system::ptr(
					new lin_constraint_system(*invariant));
	// @todo choose a domain for the solver and stick with it (use it everywhere)
	reordered_invariant->reorder(xinit.domain());
	invariant = reordered_invariant;

	for (int i = 0; i < gn; i++) {
		typename lin_constraint_system::ptr reordered_guard =
				typename lin_constraint_system::ptr(
						new lin_constraint_system(*guards[i]));
		// @todo choose a domain for the solver and stick with it (use it everywhere)
		reordered_guard->reorder(xinit.domain());
		guards[i] = reordered_guard;
	}

	z[0].begin = invariant->begin();
	z[0].end = invariant->end();
	indices[0] = 0;
	indices[1] = 2 * invariant->size();

	for (int i = 0; i < gn; i++) {
		LOGGER_OS(DEBUG7,"state_jump_intervals_solver")
			<< "guard " << i << ":\n" << guards;
		//setting z for this guard
		z[i + 1].begin = guards[i]->begin();
		z[i + 1].end = guards[i]->end();

		indices[(i + 1) * 2] = indices[(i + 1) * 2 - 1];
		indices[(i + 1) * 2 + 1] = indices[(i + 1) * 2] + 2 * guards[i]->size();
	}

	lin_constraint_state_functor z_functor(z,
			(const math::state_functor<scalar_type> *) &f, indices[2 * gn + 1]);

	/// solver initialisation
	solver.init_rootf(xinit, tinit, f, z_functor, par);
	solver.specify_jacobian(Jac);

	///do solver adjustment :
	/// a little backward then forward integration up to t = tinit in
	/// order to gather information on roots of invariants and guards at t = tinit
	/// (to detect and handle properly roots at the beginning of the algorithm)
	interm_result = solver.solve_firstrootf(tinit - 1);
	while (!approx_comparator::is_definitely_strictly_smaller(
			interm_result.stop_time, tinit)) {
		interm_result = solver.solve_firstrootf(tinit - 1);
	}
	solver.init_rootf(interm_result.stop_state, interm_result.stop_time, f,
			z_functor, par);
	while (approx_comparator::is_definitely_strictly_smaller(
			interm_result.stop_time, tinit)) {
		interm_result = solver.solve_firstrootf(tinit + par.max_timestep);
	}

	/// if solver has not stopped at t = tinit, there is no root at the beginning of the integration.
	if (approx_comparator::is_definitely_strictly_larger(
			interm_result.stop_time, tinit) && interm_result.root_info.size1()> 0) {
		interm_result.root_info = root_matrix(1, indices[gn * 2 + 1], rsNO_ROOT);
	}

	///for verfication of invariant
	/// sat is true if the invariant is true, atom is true is the invariant is atomic, e.g. only satisfied at t=tinit
	bool invariant_satisfied = true;
	bool invariant_atomic = false;

	/// root_info matrix should not be empty in the next loops
	if (interm_result.root_info.size1() == 0)
		interm_result.root_info
				= root_matrix(1, indices[gn * 2 + 1], rsNO_ROOT);

	/// evaluating invariant and then guards values and directions (if roots)  at t=tinit
	/// Goal is to detect if a guard/ the invariant is satisfied and if it is the case, verify
	/// if this situation is atomic, e.g. one of the constraint is rooted and goes up, which
	/// means that the constraint set will not be satisfied after t = current time
	/// The same sort of check is performed in the main iteration each time the solver stops with a root.
	it = invariant->begin();
	for (int i = 0; (2 * i) < indices[1]; i++) {
		LOGGER_OS(DEBUG7,"state_jump_intervals_solver")
			<< "initial check on invariant " << i << ":" << *it << " : ";
		// Check whether x violates the constraint or is on the border
		if (interm_result.root_info(0, 2 * i) == rsFROM_BELOW) {
			invariant_atomic = true;
		} else if (interm_result.root_info(0, 2 * i) == rsNO_ROOT) {
			xp = lin_constraint_evaluation(*it, xinit);
			if (approx_comparator::is_definitely_strictly_pos(xp)) {
				invariant_satisfied = false;
			} else if (!approx_comparator::is_definitely_strictly_neg(xp)) {
				// The invariant constraint is on the border.
				// Check the derivative to see the whether the invariant function is increasing or decreasing.
				if (interm_result.root_info(0, 2 * i + 1) == rsFROM_BELOW) {
					invariant_atomic = true;
				}
				/*else{
				 printf(
				 "WARNING : No root found by solver, but constraint %u is rooted",
				 i);
				 }*/
			}
		}
		++it;
	}

	/// If the invariant is not satisfied at the beginning, no trajectory should be computed and no guards taken
	/// (state is invalid !!!)
	if (!invariant_satisfied) {
		LOGGER_OS(DEBUG5,"state_jump_intervals_safe")
					<< "Invariant is not satisfied at the beginning of the simulation !!!";
		res.traj = trajectory(xinit.domain());
		res.ji = std::vector<interval_list>(gn);
		res.horizont_reached = false;
		return res;
	}

	/// If the invariant is not atomic at the beginning, no further trajectory will be computed but guards may be taken.
	/// See condition of main loop
	if (invariant_atomic)
		LOGGER_OS(DEBUG5,"state_jump_intervals_safe")
			<< "Inv atomic at beginning !!!";

	///now evaluating if guards are satisfied at beginning.
	/// the atomic et satisfied booleans play the same role as atom & satisf play for the infariant
	res.traj = trajectory(xinit.domain(), tinit, xinit.get_vector());
	res.ji = std::vector<interval_list>(gn);

	for (int i = 0; i < gn; i++) {
		it = guards[i]->begin();
		remaining_to_satisf[i] = 0;
		guard_atomic = false;
		guard_satisfied = true;
		for (int j = 0; j < guards[i]->size(); j++) {
			LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
				<< "initial check on guard " << i << ":" << *it << " : ";

			if (interm_result.root_info(0, indices[2 * (i + 1)] + 2 * j)
					== rsFROM_BELOW) {
				LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
					<< "atomic" << "\n";
				guard_atomic = true;
				remaining_to_satisf[i]++;

			} else if (interm_result.root_info(0, indices[2 * (i + 1)] + 2 * j)
					== rsNO_ROOT) {
				xp = lin_constraint_evaluation(*it, xinit);
				scalar_type inh_coeff = it->get_canonic_inh_coeff();
				if (definitely(is_GT(xp - inh_coeff, -inh_coeff))) {
					LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
						<< "violated by " << xp << "\n";

					guard_satisfied = false;
					remaining_to_satisf[i]++;
				} else if (!definitely(is_LT(xp - inh_coeff, -inh_coeff))) {
					LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
						<< "zero at " << xp << "\n";

					// Check the derivative to see the whether the guard function is increasing or decreasing.
					// The derivative is f(x). If a^Tf(x)>0, then the function is increasing.
					if (interm_result.root_info(0,
							indices[2 * (i + 1)] + 2 * j + 1) == rsFROM_BELOW) {
						LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
							<< "from below, so atomic " << "\n";
						guard_atomic = true;
						// The constraint is violated for all t'>t, so
						// it needs to be added to the constraints remaining to satisfy
						remaining_to_satisf[i]++;

					}

				}

			}

			++it;
		}
		if (guard_satisfied && guard_atomic) {
			jump_ok[i] = false;
			temp.interv.set_lower(tinit);
			temp.interv.set_upper(tinit);
			temp.index_low = 0;
			temp.index_high = 0;
			res.ji[i].push_back(temp);
		} else if (guard_satisfied) {
			jump_ok[i] = true;
			tdeb_k[i] = tinit;
			tdeb_index[i] = 0;
		} else
			jump_ok[i] = false;
	}

	solver.init_rootf(xinit, tinit, f, z_functor, par); // init solver
	solver.specify_jacobian(Jac);

	//algorithm main loop
	res.horizont_reached = false;
	bool inv = !invariant_atomic;
	scalar_type t = tinit;

	scalar_type tmax_next;
	scalar_type tprec;
	// inv && tmax>=t
	while (inv && (tmax < 0 || !maybe(approx_comparator::is_LE(tmax, t)))) {

		//TODO find a better solution if tmax < 0
		tmax_next = (tmax < 0) ? (2 * t + 1) : tmax;

		interm_result = solver.solve_firstrootf(tmax_next);
		res.traj.insert(interm_result.traj);
		tprec = t;
		t = interm_result.stop_time;

		if (!(interm_result.root_info.size1() == 0)
				&& approx_comparator::is_definitely_strictly_larger(t, tprec)) {

			LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
				<< "Root infos : NO_ROOT=" << (int) rsNO_ROOT
						<< ", FROM_ABOVE=" << (int) rsFROM_ABOVE
						<< ", FROM_BELOW=" << (int) rsFROM_BELOW << std::endl
						<< "root matrix returned by solver "
						<< interm_result.root_info << "and g values : "
						<< z_functor.map(interm_result.stop_state);

			it = invariant->begin();
			for (int i = 0; (2 * i) < indices[1]; i++) {
				xp = lin_constraint_evaluation(*it, interm_result.stop_state);
				if (interm_result.root_info(0, 2 * i) == rsFROM_BELOW) {
					inv = false;
					LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
						<< "invariant rooted, index = " << 2 * i
								<< "with value " << xp;
				} else LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
					<< "invariant not rooted, index = " << 2 * i
							<< "with value " << xp;
				++it;
			}

			for (int i = 0; i < gn; i++) //for all guards
			{

				it = guards[i]->begin();
				guard_atomic = false;
				temp_satisfied = 0;
				for (int j = 0; j < guards[i]->size(); j++) {
					LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
						<< "check on guard " << i << ", index=" << (indices[2
								* (i + 1)] + 2 * j) << " :" << *it << " : ";

					// Note: The solver seems to have quite unreliable root info
					// (timesteps too big?)
					// So let's try using the derivative

					bool below = false;
					bool above = false;

					below = interm_result.root_info(0,
							indices[2 * (i + 1)] + 2 * j) == rsFROM_BELOW;
					above = interm_result.root_info(0,
							indices[2 * (i + 1)] + 2 * j) == rsFROM_ABOVE;

					if (below) {
						LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
							<< "from below, so atomic";

						guard_atomic = true;
						remaining_to_satisf[i]++;
						temp_satisfied++;

					} else if (above) {
						remaining_to_satisf[i]--;

						LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
							<< "from above, remaining:"
									<< remaining_to_satisf[i];
					} else {

						xp = lin_constraint_evaluation(*it,
								interm_result.stop_state);
						scalar_type inh_coeff = it->get_canonic_inh_coeff();
						double dxp = lin_constraint_evaluation(*it,
								f.map(interm_result.stop_state));
						LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
							<< "solver says not root, guard value " << xp
									<< ", checking derivative dxp=" << dxp;

						if (interm_result.root_info(0,
								indices[2 * (i + 1)] + 2 * j + 1)
								== rsFROM_BELOW && is_MEQ(xp - inh_coeff,
								-inh_coeff)) {
							LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
								<< "is satisfied";
							temp_satisfied++;
						}

					}

					++it;
				}

				if (jump_ok[i] && remaining_to_satisf[i] != 0) {

					temp.interv.set_lower(tdeb_k[i]);
					temp.interv.set_upper(t);
					temp.index_low = tdeb_index[i];
					temp.index_high = res.traj.size() - 1;
					res.ji[i].push_back(temp);

					jump_ok[i] = false;

				} else if (!jump_ok[i] && remaining_to_satisf[i] == 0) {
					jump_ok[i] = true;
					tdeb_k[i] = t;
					tdeb_index[i] = res.traj.size() - 1;

				} else if (!jump_ok[i] && remaining_to_satisf[i]
						- temp_satisfied == 0) {
					temp.interv.set_lower(t);
					temp.interv.set_upper(t);
					temp.index_low = res.traj.size() - 1;
					temp.index_high = temp.index_low;
					res.ji[i].push_back(temp);
				}

			}

		}
	}

	//if(t>=tmax) res.horizont_reached=true; -->tribool
	if ((!(tmax < 0)) && maybe(approx_comparator::is_LE(tmax, t)))
		res.horizont_reached = true;

	for (int i = 0; i < gn; i++) {
		if (jump_ok[i]) {
			temp.interv.set_lower(tdeb_k[i]);
			temp.interv.set_upper(t);
			temp.index_low = tdeb_index[i];
			temp.index_high = res.traj.size() - 1;
			res.ji[i].push_back(temp);
		}
	}

	return res;

}

template<typename scalar_type>
typename traj_simu<scalar_type>::state_jump_results traj_simu<scalar_type>::state_jump_intervals_derivative(
		const state_functor & f, const state_matrix_functor & Jac,
		const state & xinit, scalar_type tinit, scalar_type tmax,
		std::vector<typename lin_constraint_system::const_ptr> & guards,
		typename lin_constraint_system::const_ptr orig_invariant) const {

	//	reorder_state_functor<scalar_type> f= reorder_state_functor<scalar_type>(f_i)
	assert(xinit.domain()==f.map(xinit).domain());

	using namespace math;
	using namespace numeric;

	unsigned int gn = guards.size();
	state_jump_results res;
	scalar_type tdeb_k[gn]; //if jump_ok[i] then tdeb_k[i] time the root occurred
	unsigned int tdeb_index[gn]; // GF: the index in the current trajectory?

	bool jump_ok[gn]; //is the jump OK at current time ?
	scalar_type xp; //scalar used for initialisation
	state xinitprim;
	scalar_type zero(0);

	ode_solver solver;

	// ---------------------------------------------------------------------
	// prepare the affine maps corresponding to the root functions
	// ---------------------------------------------------------------------
	// every guard consists of a set of constraints.
	// every constraint is considered a root function.
	// the state is in the guard if all root functions are nonpositive.
	// z_map[0] is the invariant
	// z_map[1] ... z_map[gn] are the guards
	std::vector<typename affine_map::const_ptr> z_map(gn + 1);
	// @todo choose a domain for the solver and stick with it (use it everywhere)
	// current : domain of xinit chosen.

	// to bring lin_constraint_system to domain of solver
	// + convert to affine map form for functor
	matrix<scalar_type> A;
	vector<scalar_type> b;
	positional_vdomain dom;
	positional_vdomain codom;

	// Note: The canonic matrix form is Ax<=b, to the associated root function
	// must be the affine map Ax-b==0.

	lin_constraint_system invariant(*orig_invariant);
	invariant.reorder(xinit.domain());
	canonic_matrix_form(A, b, dom, invariant);
	codom = positional_vdomain();
	for (int k = 0; k < invariant.size(); k++)
		codom.add_variable(variable("traj_simu_dummy" + to_string(k)));
	z_map[0] = typename affine_map::ptr(
			new affine_map(vdom_matrix<scalar_type> (codom, dom, A),
					state(codom, -b)));

	// Reorder the guard constraints and prepare the root functions
	for (int i = 0; i < gn; i++) {
		typename lin_constraint_system::ptr reordered_guard =
				typename lin_constraint_system::ptr(
						new lin_constraint_system(*(guards[i])));
		reordered_guard->reorder(xinit.domain());
		guards[i] = reordered_guard;

		canonic_matrix_form(A, b, dom, *(guards[i]));
		codom = positional_vdomain();
		for (int k = 0; k < guards[i]->size(); k++)
			codom.add_variable(variable("traj_simu_dummy" + to_string(k)));
		z_map[i + 1] = typename affine_map::ptr(
				new affine_map(vdom_matrix<scalar_type> (codom, dom, A),
						state(codom, -b)));
	}
	// ---------------------------------------------------------------------


	// checking invariant + which guard is satified and which is'nt.
	// guards are assumed of the form ax<=b (calcul = ax - b<=O)

	// verification of invariant
	bool invariant_satisfied = true;
	bool invariant_atomic = false;

	// Check for each constraint in the invariant whether it is satisfied or not
	for (typename lin_constraint_system::const_iterator it = invariant.begin(); invariant_satisfied
			&& it != invariant.end(); ++it) {
		// Check whether x violates the constraint or is on the border
		xp = it->normal_eval(xinit);
		// For better numerics, don't compare xp to zero, but
		// xp-inh_coeff to -inh_coeff. That way we take into account the relative error
		scalar_type inh_coeff = it->get_canonic_inh_coeff();
		if (definitely(is_GT(xp, -inh_coeff))) {
			// The invariant constraint is violated.
			LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
				<< "inv violated with " << xp << " > " << -inh_coeff;
			invariant_satisfied = false;
		} else if (!definitely(is_LT(xp, -inh_coeff))) {
			// The invariant constraint is on the border.
			// Check the derivative to see the whether the invariant function is increasing or decreasing.
			// The derivative is f(x). If a^Tf(x)>0, then the function is increasing, and so the invariant is atomic.
			scalar_type dxp = it->normal_eval(f.map(xinit));
			if (definitely(is_GT(dxp, scalar_type(0)))) {
				LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
					<< "inv from below, so atomic " << "\n";
				invariant_atomic = true;
			}
			// Note: for now we don't check whether the derivative is zero,
			// which could cause problems with the solver
		}
	}
	// If the invariant is violated from the start, return an empty trajectory
	if (!invariant_satisfied) {
		LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
			<< ("Inv false at beginning !!!\n");
		res.traj = trajectory(xinit.domain());
		res.ji = std::vector<interval_list>(gn);
		res.horizont_reached = false;
		return res;
	}

	if (invariant_atomic)
		LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
			<< ("Inv atomic at beginning !!!\n");

	// We found our first point, let's add it to the result trajectory
	res.traj = trajectory(xinit.domain(), tinit, xinit.get_vector());
	res.ji = std::vector<interval_list>(gn);

	// currently sought roots
	// if z[i] < 0 then all constraints of inv/guard need to be evaluated
	// otherwise only the constraint with index z[i]
	std::vector<int> z(gn + 1);
	// for the guard i, the root function values of the constraints in the
	// interval [index_beg[i+1],index_end[i+1]-1] need to be verified
	// These indices correspond to the position in the output vector of
	// z_functor.
	int index_beg[gn + 1];
	int index_end[gn + 1];
	// all constraints of the invariant need to be tracked
	z[0] = -1;
	index_beg[0] = 0;
	index_end[0] = invariant.size();

	// analyze all of the guards on their status, and find which
	// constraints need to be evaluated
	for (int i = 0; i < gn; i++) {
		scalar_type max_val(0);
		bool guard_satisfied = true;
		bool guard_atomic = false;
		int position = 0;
		int pos_max = -1;
		int pos_atom = -1;

		// find out if the guard is satisfied, atomic, and the "farthest" constraint
		for (typename lin_constraint_system::const_iterator it =
				guards[i]->begin(); it != guards[i]->end(); ++it) {
			scalar_type xp = it->normal_eval(xinit);
			scalar_type inh_coeff = it->get_canonic_inh_coeff();
			if (pos_max < 0 || max_val < xp + inh_coeff) {
				pos_max = position;
				max_val = xp + inh_coeff;
			}

			if (guard_satisfied) {
				if (definitely(is_GT(xp, -inh_coeff))) {
					LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
						<< "guard " << i << " constraint " << *it << " violated by " << xp + inh_coeff << "\n";
					guard_satisfied = false;
				} else if (!definitely(is_LT(xp, -inh_coeff))) {
					LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
						<< "guard " << i << " constraint " << *it << " zero at " << xp + inh_coeff << "\n";
					// Check the derivative to see the whether the guard function is increasing or decreasing.
					// The derivative is f(x). If a^Tf(x)>0, then the function is increasing.
					scalar_type dxp = it->normal_eval(f.map(xinit));
					//				if(approx_comparator::is_maybe_equal(xp-inh_coeff,-inh_coeff) )
					//					std::runtime_error("Root found with rooted derivative, don't know how to deal with !");
					if (definitely(is_GT(dxp, scalar_type(0)))) {
						LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
							<< "from below, so atomic " << "\n";
						guard_atomic = true;

						pos_atom = position;
					}
				}
			}

			position++;
		}

		if (guard_satisfied && guard_atomic) {
			jump_interval temp;
			temp.interv.set_lower(tinit);
			temp.interv.set_upper(tinit);
			temp.index_low = 0;
			temp.index_high = 0;
			res.ji[i].push_back(temp);

			z[i + 1] = pos_atom;
			jump_ok[i] = false;

			// register the position in the output vector
			index_beg[i + 1] = index_end[i]; // continue after the index of the previous constraint
			index_end[i + 1] = index_beg[i + 1] + 1;
		} else if (guard_satisfied) {

			z[i + 1] = -1;
			jump_ok[i] = true;

			tdeb_k[i] = tinit;
			tdeb_index[i] = 0;

			// register the position in the output vector
			index_beg[i + 1] = index_end[i];
			index_end[i + 1] = index_beg[i + 1] + guards[i]->size(); // reserve room for all constraints
		} else {
			LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
				<< "guard " << i << " violated, adding constraint " << pos_max << "\n";
			z[i + 1] = pos_max;
			jump_ok[i] = false;

			// register the position in the output vector
			index_beg[i + 1] = index_end[i]; // continue after the index of the previous constraint
			index_end[i + 1] = index_beg[i + 1] + 1;
		}

	}

	//lin_constraint_state_functor z_functor(z,(const math::state_functor<scalar_type> *)&f,indices[(gn<<1)+1]);
	affine_map_vector_state_vector_functor z_functor(z_map, z, index_end[gn]);
	solver.init_rootf(xinit, tinit, f, z_functor, par); // init solver
	solver.specify_jacobian(Jac);

	// end of init

	typename ode_solver::rootf_result interm_result;

	typedef typename ode_solver::rootf_result::root_status root_status;

	root_status rsNO_ROOT = ode_solver::rootf_result::NO_ROOT;
	root_status rsFROM_BELOW = ode_solver::rootf_result::FROM_BELOW;
	root_status rsFROM_ABOVE = ode_solver::rootf_result::FROM_ABOVE;

	//algorithm body
	res.horizont_reached = false;
	bool inv = !invariant_atomic;
	scalar_type t = tinit;
	bool at_least_one_from_below;
	bool solver_crash = false;

	scalar_type tmax_next;
	// inv && tmax>=t
	//	while(inv && ( tmax<0 || !maybe(approx_comparator::is_LE(tmax,t) ) ) && !solver_crash )
	while (inv && (tmax < 0 || definitely(approx_comparator::is_LE(t, tmax)))
			&& !solver_crash) {
		tmax_next = (tmax < 0) ? (2 * t + 1) : tmax;

		z_functor.update_approx_size(index_end[gn]);
		LOGGER_OS(DEBUG7,"state_jump_intervals_derivative")
				<< "running solver with "<< "guards: " << z_map << "roots: " << z;
		try {
			solver.change_rootf(z_functor);
			interm_result = solver.solve_firstrootf(tmax_next);
			res.traj.insert(interm_result.traj);
			t = interm_result.stop_time;
		} catch (std::exception& e) {
			basic_warning("state_jump_intervals_derivative",
					"ODE solver crashed, result is incomplete.",
					basic_warning::INCOMPLETE_OUTPUT);
			solver_crash = true;
			IFLOGGER(DEBUG4) {
				LOGGER_OS(DEBUG4,"state_jump_intervals_derivative")
						<< "guards: " << z_map << "roots: " << z;
				throw e;
			}
		}

		if (!(interm_result.root_info.size1() == 0)) {

			LOGGER_OS(DEBUG7,"state_jump_intervals_derivative")
				<< " stop_state : " << interm_result.stop_state
						<< "and inv & g values : " << z_functor.map(
						interm_result.stop_state);

			{
				typename lin_constraint_system::const_iterator it =
						invariant.begin();
				for (int i = index_beg[0]; i < index_end[0]; i++) {
					xp = it->normal_eval(interm_result.stop_state);
					scalar_type inh_coeff = it->get_canonic_inh_coeff();
					scalar_type dxp = it->normal_eval(
							f.map(interm_result.stop_state));
					//				if(approx_comparator::is_maybe_equal(xp-inh_coeff,-inh_coeff) && approx_comparator::is_maybe_equal(dxp-inh_coeff,-inh_coeff) )
					//					std::runtime_error("Root found with rooted derivative, don't know how to deal with !");

					// The state is outside the invariant, or on the border of the invariant with derivative pointing outside
					if (definitely(is_GT(xp, -inh_coeff)) || (is_MEQ(xp,
							-inh_coeff) && definitely(
							is_GT(dxp, scalar_type(0))))) {
						inv = false;
					}
					++it;
				}
			}

			for (int i = 0; i < gn; i++) //for all guards
			{
				typename lin_constraint_system::const_iterator it =
						guards[i]->begin();
				if (jump_ok[i]) {
					// the state was previously inside the guard
					int j, git, posgit;
					git = 0; // index counter
					posgit = -1; // index of the positive constraint
					for (j = index_beg[i + 1]; jump_ok[i] && j < index_end[i
							+ 1]; ++j) {
						// here we could use j to access the output vector of z_functor.
						// for now, we re-evaluate the constraints
						xp = it->normal_eval(interm_result.stop_state);
						scalar_type inh_coeff = it->get_canonic_inh_coeff();
						scalar_type dxp = it->normal_eval(
								f.map(interm_result.stop_state));
						//						if(approx_comparator::is_maybe_equal(xp-inh_coeff,-inh_coeff) && approx_comparator::is_maybe_equal(dxp-inh_coeff,-inh_coeff) )
						//							std::runtime_error("Root found with rooted derivative, don't know how to deal with !");
						// the state is on the border of the guard, and the derivative is pointing outside
						if (definitely(is_GT(xp, -inh_coeff)) || (is_MEQ(xp,
								-inh_coeff) && definitely(
								is_GT(dxp, scalar_type(0))))) {
							jump_ok[i] = false;
							posgit = git;
						}

						++git;
						++it;
					}
					if (jump_ok[i] == false) // change tested constraint if violation has been found
					{
						z[i + 1] = posgit;
						index_beg[i + 1] = index_end[i];
						index_end[i + 1] = index_beg[i + 1] + 1;

						jump_interval temp;
						temp.interv.set_lower(tdeb_k[i]);
						temp.interv.set_upper(t);
						temp.index_low = tdeb_index[i];
						temp.index_high = res.traj.size() - 1;
						res.ji[i].push_back(temp);
					} else // tested constraints remain the same, just updating indexes
					{
						index_beg[i + 1] = index_end[i];
						index_end[i + 1] = index_beg[i + 1] + guards[i]->size();
					}
				} else {
					// the state was outside of this guard
					// check if the constraint is still violated

					// move iterator to the constraint with index z[i+1]
					typename lin_constraint_system::const_iterator it =
							guards[i]->begin();
					for (int kk = 0; kk < z[i + 1]; ++kk) {
						++it;
					}
					scalar_type xp = it->normal_eval(interm_result.stop_state);
					scalar_type inh_coeff = it->get_canonic_inh_coeff();

					// if not violated, restart
					if (maybe(is_LE(xp, -inh_coeff))) {
						scalar_type max_val(0);
						bool guard_satisfied = true;
						bool guard_atomic = false;
						int position = 0;
						int pos_max = -1;
						int pos_atom = -1;

						// find out if the guard is satisfied, atomic, and the "farthest" constraint
						for (typename lin_constraint_system::const_iterator it =
								guards[i]->begin(); it != guards[i]->end(); ++it) {
							scalar_type xp = it->normal_eval(
									interm_result.stop_state);
							scalar_type inh_coeff = it->get_canonic_inh_coeff();
							if (pos_max < 0 || max_val < xp + inh_coeff) {
								pos_max = position;
								max_val = xp + inh_coeff;
							}

							if (guard_satisfied) {
								if (definitely(is_GT(xp, -inh_coeff))) {
									LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
										<< "violated by " << xp + inh_coeff
												<< "\n";
									guard_satisfied = false;
								} else if (!definitely(is_LT(xp, -inh_coeff))) {
									LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
										<< "zero at " << xp + inh_coeff << "\n";
									// Check the derivative to see the whether the guard function is increasing or decreasing.
									// The derivative is f(x). If a^Tf(x)>0, then the function is increasing.
									scalar_type dxp = it->normal_eval(
											f.map(interm_result.stop_state));
									//				if(approx_comparator::is_maybe_equal(xp-inh_coeff,-inh_coeff) )
									//					std::runtime_error("Root found with rooted derivative, don't know how to deal with !");
									if (definitely(is_GT(dxp, scalar_type(0)))) {
										LOGGER_OS(DEBUG7,"state_jump_intervals_safe")
											<< "from below, so atomic " << "\n";
										guard_atomic = true;

										pos_atom = position;
									}
								}
							}

							position++;
						}
						// Note: guard_atomic may be true if guard_satisfied is false, because
						// guard_satisfied can change after guard_atomic has become true

						if (guard_satisfied && guard_atomic) {
							jump_interval temp;
							temp.interv.set_lower(t);
							temp.interv.set_upper(t);
							temp.index_low = res.traj.size() - 1;
							temp.index_high = temp.index_low;
							res.ji[i].push_back(temp);

							z[i + 1] = pos_atom;
							jump_ok[i] = false;

							// register the position in the output vector
							index_beg[i + 1] = index_end[i]; // continue after the index of the previous constraint
							index_end[i + 1] = index_beg[i + 1] + 1;
						} else if (guard_satisfied) {

							z[i + 1] = -1;
							jump_ok[i] = true;

							tdeb_k[i] = t;
							tdeb_index[i] = res.traj.size() - 1;

							// register the position in the output vector
							index_beg[i + 1] = index_end[i];
							index_end[i + 1] = index_beg[i + 1]
									+ guards[i]->size(); // reserve room for all constraints
						} else {
							z[i + 1] = pos_max;
							jump_ok[i] = false;

							// register the position in the output vector
							index_beg[i + 1] = index_end[i]; // continue after the index of the previous constraint
							index_end[i + 1] = index_beg[i + 1] + 1;
						}
					} else // tested constraints remain the same, just updating indexes
					{
						index_beg[i + 1] = index_end[i];
						index_end[i + 1] = index_beg[i + 1] + guards[i]->size();
					}
				} // else: state was outside of guard

			} // for all guards

		} // if there is a root

	}

	//if(t>=tmax) res.horizont_reached=true; -->tribool
	if ((!(tmax < 0)) && maybe(is_LE(tmax, t)))
		res.horizont_reached = true;

	for (int i = 0; i < gn; i++) {
		if (jump_ok[i]) {
			jump_interval temp;
			temp.interv.set_lower(tdeb_k[i]);
			temp.interv.set_upper(t);
			temp.index_low = tdeb_index[i];
			temp.index_high = res.traj.size() - 1;
			res.ji[i].push_back(temp);
		}
	}

	return res;

}

template<typename scalar_type>
typename traj_simu<scalar_type>::roots_and_trajectory traj_simu<scalar_type>::select_next_traj_states(
		const state_jump_results & res) const {

	unsigned int mp;
	roots_and_trajectory ret;
	if (res.horizont_reached) {
		static bool has_warned(false);
		if (has_warned == false) {
			has_warned = true;
			basic_warning("traj_simu<scalar_type>::select_next_traj_states",
					"Reached time horizon without exhausting all states, result is incomplete.",
					basic_warning::INCOMPLETE_OUTPUT);
		}
	}

	//TODO improve this code ?
	interval_list l;
	ret.roots = std::vector<trajectory>();
	ret.roots.reserve(res.ji.size());

	for (int i = 0; i < res.ji.size(); i++) {
		l = res.ji[i];
		trajectory ts;
		//ts=trajectory(res.traj.domain());
		for (typename interval_list::const_iterator it = l.begin(); it
				!= l.end(); ++it) {
			ts.insert((*it).interv.lower().get_val(),
					res.traj.get_state((*it).index_low));
			if ((*it).index_high != (*it).index_low)
				ts.insert((*it).interv.upper().get_val(),
						res.traj.get_state((*it).index_high).get_vector());
			/* GF: Taking out the middle point; it doesn't seem to
			 * be very useful and increases complexity.
			 */
			/*
			mp = ((*it).index_high - (*it).index_low) / 2;
			if (mp > 0) {
				mp = mp + (*it).index_low;
				ts.insert(res.traj.get_time(mp), res.traj.get_state(mp));
			}
			*/
		}
		LOGGER_OS(DEBUG7,"select_next_traj_states")
			<< "selected roots:\n" << ts;
		ret.roots.push_back(ts);
	}
	ret.traj = trajectory(res.traj);
	assert ( ret.roots.size()==res.ji.size() );
	return ret;

}

template<typename scalar_type>
typename traj_simu<scalar_type>::simulation_results traj_simu<scalar_type>::image_next_traj_states(
		const simulation_results & res,
		const std::vector<typename state_functor::const_ptr> & tr_assign) const {

	assert(res.size()==tr_assign.size());

	//this code is a little bit messy, and creates a lot of new objects.
	// perhaps use a visitor thing  to avoid this?

	simulation_results ret = simulation_results();
	int j = 0;
	for (typename simulation_results::const_iterator it = res.begin(); it
			!= res.end(); ++it) {

		if (it->size() > 0) {
			// Note: the domain after a jump can be different than before.
			trajectory t;
			const math::vector<scalar_type>& times = it->get_times();
			//		reorder_state_functor<scalar_type> rftemp = reorder_state_functor<scalar_type>(*tr_assign[j]);
			for (int i = 0; i < it->size(); i++) {
				t.insert(
						times[i],
						state_image_through_transition(it->get_state(i),
								*tr_assign[j]));
			}
			LOGGER_OS(DEBUG7,"image_next_traj_states")
				<< "mapped roots:\n" << t;
			ret.push_back(t);
		} else
			ret.push_back(trajectory());
		++j;
	}

	return ret;
}

template<typename scalar_type>
typename traj_simu<scalar_type>::state traj_simu<scalar_type>::state_image_through_transition(
		const state & s, const state_functor & f_assign) const {

	return f_assign.map(s);
}

#endif
