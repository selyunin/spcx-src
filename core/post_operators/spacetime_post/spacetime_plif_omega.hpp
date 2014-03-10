/*
 * spacetime_plif_omega.hpp
 *
 *  Created on: Oct 20, 2012
 *      Author: notroot
 */

#ifndef SPACETIME_PLIF_OMEGA_HPP_
#define SPACETIME_PLIF_OMEGA_HPP_

#include "spacetime_plif.h"

#include <stack>
#include "utility/basic_exception.h"

namespace spacetime {

inline
spacetime_plif::timed_directional_interval_sequence spacetime_plif::omega_X0_err(
		const direction& d, const duration& t, const duration& delta,
		const error_type& err_begin, const error_type& err_end) {
	LOGGERSW(DEBUG6,"omega_X0_err","computing evolution of autonomous dynamics");

//	std::cout << "omega_X0 in direction " << d << " delta=" << delta;
	//			<< " with A=" << my_aff.dyn.get_A() << ", exp(Adelta)="
	//			<< exp_AdeltaT(delta).transpose() << std::endl;

	//	vector_type d_end = exp_AdeltaT(delta) * d;
	vector_type vec_Omega_beg = rho_vec_X0(d);
	//	vector_type vec_Omega_end = rho_vec_X0(d_end);
	scalar_type rho_Omega_beg = scalar_product(d, vec_Omega_beg);
	//	scalar_type rho_Omega_end = scalar_product(d_end, vec_Omega_end);

	// get the sequence of error terms
	vector_type vec_Omega_end;  // unused return argument
	timed_directional_interval_sequence seq = omega_X0_delta(d, t, delta,
			err_begin, err_end, vec_Omega_beg, vec_Omega_end);

//	std::cout << "Omega_X0_err sequence: " << std::endl;
//	seq.print(std::cout);

	// add initial element
	timed_directional_interval tdi;
	tdi.t = t;
	tdi.d = d;
	if (seq.begin() != seq.end()) {
		tdi.delta = seq.begin()->t - tdi.t;
	} else {
		tdi.delta = delta;
	}
	tdi.itv = scalar_interval(rho_Omega_beg,rho_Omega_beg);
	seq.push_front(tdi);

	// add last element
	//	tdi.t = t + delta;
	//	tdi.d = d_end;
	//	tdi.delta = duration(0);
	//	tdi.itv = scalar_interval(rho_Omega_end);
	//	seq.push_back(tdi);

//	std::cout << "Omega_X0_delta sequence: " << std::endl;
//	seq.print(std::cout);

	return seq;
}

inline
spacetime_plif::timed_directional_interval_sequence spacetime_plif::omega_X0_err_simplified(
		const direction& d_ini, const duration& t_ini, const duration& delta_ini,
		const error_type& err_begin_ini, const error_type& err_end_ini) {
	LOGGERSW(DEBUG6,"omega_X0_err","computing evolution of autonomous dynamics");

	using namespace math;

	timed_directional_interval tdi;
	timed_directional_interval_sequence seq;

	direction d = d_ini;
	duration t = t_ini;
	duration delta = delta_ini;

	error_type err_begin = err_begin_ini;
	error_type err_end = err_end_ini;
	duration t_max = t_ini + delta_ini;

	vector_type d_end, ef, eb;
	vector_type vec_Omega_beg, vec_Omega_end;
	scalar_type rho_Omega_beg, rho_Omega_end;

	// Do the initial computation
	//	vector_type d_end = exp_AdeltaT(delta) * d;
	vec_Omega_beg = rho_vec_X0(d);
	//	vector_type vec_Omega_end = rho_vec_X0(d_end);
	rho_Omega_beg = scalar_product(d, vec_Omega_beg);

	// Use a LIFO queue
	std::stack<omega_X0_err_state> state_queue;
	omega_X0_err_state queue_state;

	bool err_too_big;
	bool skip_to_smaller_delta;

	scalar_interval itv_beg,itv_end;

	scalar_type lower_bound,upper_bound;
	bool use_max_slope_forw;
	scalar_type max_slope_forw;
	bool use_min_slope_back;
	scalar_type min_slope_back;
	scalar_type rho_slope;

	do {
		// Find a suitable delta_curr that is low enough numerically and meets error bounds
		do {
			err_too_big = false;
			skip_to_smaller_delta = (min_infeasible_delta > scalar_type(0)
					&& delta >= min_infeasible_delta / scalar_type(16));

			if (!skip_to_smaller_delta) {
				// try to obtain the parameters for *this; if there are numerical issues
				// skip to next recursion with smaller delta
				try {
					d_end = exp_AdeltaT(delta) * d;
					// get the vectors for E_Omega(delta,t)
					ef = math::vec_abs(rho_vec(d, E_X0f(delta)));
					eb = math::vec_abs(rho_vec(d, E_X0b(delta)));
				} catch (math::invalid_number_exception& e) {
					LOGGER(DEBUG6, __FUNCTION__,
							"decreasing current time step " + to_string(delta)
									+ " because of numerical issues");
					skip_to_smaller_delta = true;
				}
//	std::cout << std::endl	<< "taking delta " << delta << std::endl;
			}
			if (!skip_to_smaller_delta) {
				// compute the support
				vec_Omega_end = rho_vec_X0(d_end);
				rho_Omega_end = scalar_product(d_end, vec_Omega_end);

				scalar_type l_absc_forw = rho_Omega_beg;
				scalar_type l_absc_backw = scalar_product(d, vec_Omega_end);
				scalar_type l_slope_forw = scalar_product(d_end, vec_Omega_beg)
						- l_absc_forw;
				scalar_type l_slope_backw = rho_Omega_end - l_absc_backw;

				// track the max forward slope; init with linear interp since it's known to bound from below
				// heuristic: use only at start
				use_max_slope_forw = is_MEQ(t,t_ini);
				use_min_slope_back = maybe(is_GE(t+delta,t_max));
				rho_slope = (rho_Omega_end-rho_Omega_beg)/delta;
				max_slope_forw = rho_slope;
				min_slope_back = rho_slope;

				//std::cout << "d_beg:"<<d_beg <<" d_end:" << d_end << " vec_beg:"<<  vec_Omega_beg << " vec_end:"<<  vec_Omega_end ;
//		std::cout << std::endl << "ef:" << ef << " eb:" << eb << std::endl;

				// compute the breakpoints from lambda = 0 to 1
				scalar_type lambda(0);
				scalar_type err_upp(0), err_low(0);
				scalar_type err_upp_max(0), err_low_min(0);

				// Note: deltas will be fixed later
				// tdi.delta = duration(0);
				// Note: no need to compute d, add only to first element

				// compute crossing point
				scalar_type lambda_sw = scalar_type(2.0); // start with invalid;
				if (!is_MEQ(l_slope_forw, l_slope_backw)) {
					lambda_sw = (l_absc_backw - l_absc_forw)
							/ (l_slope_forw - l_slope_backw);
					// std::cout << "lambda at switch:" << lambda << std::endl;
					if (maybe(is_LE(lambda_sw, scalar_type(0)))
							|| maybe(is_GE(lambda_sw, scalar_type(1)))) {
						lambda_sw = scalar_type(2.0); // invalid
					}
				}
				//	std::cout << "lambdasw=" << lambda_sw << " af "
				//							<< l_absc_forw << " ab " << l_absc_backw << " sl_f "
				//							<< l_slope_forw << " sl_b " << l_slope_backw << std::endl;

				// create the absolute values of the direction
				vector_type::vector_type d_abs(math::vec_abs(d.get_vector()));

				// cycle through switching lambdas
				for (size_t i = 0; i < ef.size(); ++i) {
					if (math::definitely(
							math::numeric::is_GT(ef[i] + eb[i],
									scalar_type(0)))) {
						/** Careful, lambda could be very small if the error is large */
						lambda = eb[i] / (ef[i] + eb[i]);
						// std::cout << "lambda at switch:" << lambda << std::endl;

						//err_upp = rho_E_Omega(ef, eb, d, lambda);
						// computes scalar_product(math::min(t * ef, (delta - t) * eb),  math::vec_abs(l));
						scalar_type err_upp = scalar_type(0);
						duration delta_minus_t = duration(1) - lambda;
						for (size_t i = 0; i < ef.size(); ++i) {
							scalar_type t_ef = lambda * ef[i];
							scalar_type one_t_eb = delta_minus_t * eb[i];
							err_upp += std::min(lambda * ef[i],
									delta_minus_t * eb[i]) * d_abs[i];
						}
						//std::cout << "err_upp:" << err_upp << std::endl;
						err_low = -err_upp;

						upper_bound = (scalar_type(1) - lambda)
								* rho_Omega_beg + lambda * rho_Omega_end + err_upp;
						lower_bound = std::max(
								l_absc_forw + lambda * l_slope_forw,
								l_absc_backw + lambda * l_slope_backw) + err_low;
						err_upp_max = std::max(err_upp_max, err_upp);
						// measure the lower bound error
						scalar_type err_low_meas = (scalar_type(1) - lambda)
										* rho_Omega_beg + lambda * rho_Omega_end - lower_bound;
						err_low_min = std::min(err_low_min, err_low_meas);

						// track the max slope
						if (definitely(is_GT(lambda,scalar_type(0)))) {
							max_slope_forw = std::max(max_slope_forw,(upper_bound-rho_Omega_beg)/(lambda*delta));
						} else {
							use_max_slope_forw = false;
						}
						if (definitely(is_LT(lambda,scalar_type(1)))) {
							min_slope_back = std::min(min_slope_back,(rho_Omega_end-upper_bound)/(delta-lambda*delta));
						} else {
							use_min_slope_back = false;
						}
					}
				}
				if (lambda_sw < scalar_type(1)) {
					//		std::cout << "inserting lambda_sw="<< lambda_sw << " after " << lambda << std::endl;

					err_upp = rho_E_Omega(ef, eb, d, lambda_sw);
					err_low = -err_upp;

					err_upp_max = std::max(err_upp_max, err_upp);
					err_low_min = std::min(err_low_min, err_low);

//					tdi.t = t + lambda_sw * delta;
					upper_bound = (scalar_type(1) - lambda_sw)
							* rho_Omega_beg + lambda_sw * rho_Omega_end
							+ err_upp;
					lower_bound = std::max(
							l_absc_forw + lambda_sw * l_slope_forw,
							l_absc_backw + lambda_sw * l_slope_backw)
							+ err_low;

					err_upp_max = std::max(err_upp_max, err_upp);
					// measure the lower bound error
					scalar_type err_low_meas = (scalar_type(1) - lambda_sw)
									* rho_Omega_beg + lambda_sw * rho_Omega_end - lower_bound;
					err_low_min = std::min(err_low_min, err_low_meas);
				}

				// check whether the error is within the allowed bounds:
				// err(tau) <= err_begin + (tau-t)/delta*(err_end-err_begin)

				upper_bound = rho_Omega_beg + err_upp_max;
				lower_bound = rho_Omega_beg + err_low_min;
				itv_beg = scalar_interval(lower_bound, upper_bound);
				upper_bound = rho_Omega_end + err_upp_max;
				lower_bound = rho_Omega_end + err_low_min;
				itv_end = scalar_interval(lower_bound, upper_bound);

				// compute the error
				scalar_type err = err_upp - err_low;

				err_too_big = err_too_big
						|| definitely(!err_begin.is_satisfied(itv_beg));

				err_too_big = err_too_big
						|| definitely(!err_end.is_satisfied(itv_end));
			}

			// if error too big, split
			if (err_too_big || skip_to_smaller_delta) {
				// split the interval in two parts, first and second
				// choose the time point at which to split: t+delta_first
				duration delta_split = delta / duration(2.0); // this can be any value strictly between 0 and delta
				// the allowable error bound at the split
				error_type err_split = err_begin
						+ delta_split / delta * (err_end - err_begin);

				// push the current end to the queue and compute a new one
				queue_state.t = t + delta_split;
				queue_state.delta = delta - delta_split;

				if (!skip_to_smaller_delta) {
					queue_state.d_end = d_end;
					queue_state.vec_end = vec_Omega_end;
					queue_state.rho_end = rho_Omega_end;
				}
				queue_state.err_end = err_end;
				state_queue.push(queue_state);

				// define the new time interval
				delta = delta_split;
				err_end = err_split;
			}
		} while (err_too_big || skip_to_smaller_delta);

		// accept element -- produce left-closed interval
		use_max_slope_forw = use_max_slope_forw && definitely(is_GT(max_slope_forw,rho_slope));
		use_min_slope_back = use_min_slope_back && definitely(is_LT(min_slope_back,rho_slope));
		if (!use_max_slope_forw) {
			tdi.t = t;
//		tdi.d = d_end;
			tdi.delta = delta;
			tdi.itv = itv_beg;
		} else {
			/** find the intersection of the max forward slop with
			 * the linear interpolation:
			 * rho_beg + slope*delta_inters = r1 + (r2-r1)/delta*delta_inters
			 * <=>
			 * delta_inters = (rho_beg - r1)/((r2-r1)/delta-slope)
			 */
			scalar_type delta_inters = (rho_Omega_beg - itv_beg.upper())/(rho_slope-max_slope_forw);

			// First point is a rho_beg
			tdi.t = t;
//		tdi.d = d_end;
			tdi.delta = delta_inters;
			tdi.itv = scalar_interval(rho_Omega_beg, rho_Omega_beg);
			seq.push_back(tdi);

			// Second point is at intersection
			tdi.t = t + delta_inters;
//		tdi.d = d_end;
			tdi.delta = delta - delta_inters;
			tdi.itv = itv_beg+rho_slope*delta_inters;
		}
		if (!use_min_slope_back) {
			seq.push_back(tdi);
			tdi.t = t + delta;
//		tdi.d = d_end;
			tdi.delta = duration(0);
			tdi.itv = itv_end;
			seq.push_back(tdi);
		} else {
			/** find the intersection of the min backward slope with
			 * the linear interpolation:
			 * rho_end - slope*(delta-delta_inters) = r2 - (r2-r1)/delta*(delta-delta_inters)
			 * <=>
			 * delta_inters = delta+(rho_end - r2)/((r2-r1)/delta-slope)
			 */
			scalar_type delta_inters = delta + (rho_Omega_end - itv_end.upper())/(rho_slope-min_slope_back);

			// adjust middle point before adding it
			if (definitely(is_LT(delta_inters, tdi.delta))) {
				tdi.delta -= delta_inters;
				seq.push_back(tdi);

				// add intersection of backwards slope
				tdi.t = tdi.t + tdi.delta;
//		tdi.d = d_end;
				tdi.delta = delta - tdi.t;
				tdi.itv = itv_end - rho_slope * delta_inters;
				seq.push_back(tdi);
			} else {
				// we could compute the intersection point
				// of forward and backward, but what the heck
				/** @todo Check why this is so -- it should not be possible given it's a concave function */

				// let's just keep the forward point
				seq.push_back(tdi);
			}

			// add final point
			tdi.t = t + delta;
//		tdi.d = d_end;
			tdi.delta = duration(0);
			tdi.itv = scalar_interval(rho_Omega_end, rho_Omega_end);
			seq.push_back(tdi);
		}


		// Skip to the next element
		// i.e., the end becomes the beginning of the next interval
		vec_Omega_beg = vec_Omega_end;
		rho_Omega_beg = rho_Omega_end;
		err_begin = err_end;
		d = d_end;
		t = t + delta;

		// get the next state from the queue
		if (definitely(is_LT(t, t_max))) {
			if (!state_queue.empty()) {
				queue_state = state_queue.top();
				state_queue.pop();
				t = queue_state.t;
				delta = queue_state.delta;
				d_end = queue_state.d_end;
				vec_Omega_end = queue_state.vec_end;
				rho_Omega_end = queue_state.rho_end;
				err_end = queue_state.err_end;
			} else {
				std::cerr << "t: " << t << ", t_max: " << t_max << std::endl;
				throw std::runtime_error("empty queue");
			}
		}

	} while (definitely(is_LT(t, t_max)));

	return seq;
}

inline
spacetime_plif::timed_directional_interval_sequence spacetime_plif::omega_X0_delta(
		const direction& d, const duration& t, const duration& delta,
		const error_type& err_begin, const error_type& err_end,
		const vector_type& vec_Omega_beg, vector_type& vec_Omega_end) {

	timed_directional_interval tdi;
	timed_directional_interval_sequence seq;

	const direction& d_beg = d;
	vector_type d_end, ef, eb;
	scalar_type rho_Omega_end;

	bool err_too_big = false;

	bool skip_to_smaller_delta = (min_infeasible_delta > scalar_type(0) && delta >= min_infeasible_delta/scalar_type(16));

	if (!skip_to_smaller_delta) {
	// try to obtain the parameters for *this; if there are numerical issues
	// skip to next recursion with smaller delta
	try {
		d_end = exp_AdeltaT(delta) * d;
		// get the vectors for E_Omega(delta,t)
		ef = math::vec_abs(rho_vec(d, E_X0f(delta)));
		eb = math::vec_abs(rho_vec(d, E_X0b(delta)));
	} catch (math::invalid_number_exception& e) {
		LOGGER(DEBUG6, __FUNCTION__,
				"decreasing current time step "+to_string(delta)+" because of numerical issues");
		skip_to_smaller_delta = true;
	}
//	std::cout << std::endl	<< "taking delta " << delta << std::endl;
	}
	if (!skip_to_smaller_delta) {
		vec_Omega_end = rho_vec_X0(d_end);

		scalar_type rho_Omega_beg = scalar_product(d_beg, vec_Omega_beg);
		rho_Omega_end = scalar_product(d_end, vec_Omega_end);
		scalar_type l_absc_forw = rho_Omega_beg;
		scalar_type l_absc_backw = scalar_product(d_beg, vec_Omega_end);
		scalar_type l_slope_forw = scalar_product(d_end, vec_Omega_beg)
				- l_absc_forw;
		scalar_type l_slope_backw = rho_Omega_end - l_absc_backw;

		//std::cout << "d_beg:"<<d_beg <<" d_end:" << d_end << " vec_beg:"<<  vec_Omega_beg << " vec_end:"<<  vec_Omega_end ;
//		std::cout << std::endl << "ef:" << ef << " eb:" << eb << std::endl;

		// compute the breakpoints from lambda = 0 to 1
		scalar_type lambda(0);
		scalar_type err_upp(0), err_low(0);
		scalar_type err_upp_max(0), err_low_max(0);

		// Note: deltas will be fixed later
		// tdi.delta = duration(0);
		// Note: no need to compute d, add only to first element

		// compute crossing point
		scalar_type lambda_sw = scalar_type(2.0); // start with invalid;
		if (!is_MEQ(l_slope_forw, l_slope_backw)) {
			lambda_sw = (l_absc_backw - l_absc_forw) / (l_slope_forw
					- l_slope_backw);
			// std::cout << "lambda at switch:" << lambda << std::endl;
			if (maybe(is_LE(lambda_sw, scalar_type(0))) || maybe(
					is_GE(lambda_sw, scalar_type(1)))) {
				lambda_sw = scalar_type(2.0); // invalid
			}
		}
		//	std::cout << "lambdasw=" << lambda_sw << " af "
		//							<< l_absc_forw << " ab " << l_absc_backw << " sl_f "
		//							<< l_slope_forw << " sl_b " << l_slope_backw << std::endl;

		// remember where we started inserting
		timed_directional_interval_sequence::iterator before_seq=seq.begin();

		// create the absolute values of the direction
		vector_type::vector_type d_abs(math::vec_abs(d.get_vector()));

		// cycle through switching lambdas
		for (size_t i = 0; i < ef.size(); ++i) {
			if (math::definitely(
					math::numeric::is_GT(ef[i] + eb[i], scalar_type(0)))) {
				/** Careful, lambda could be very small if the error is large */
				lambda = eb[i] / (ef[i] + eb[i]);
				// std::cout << "lambda at switch:" << lambda << std::endl;
//				if (definitely(is_GT(lambda, scalar_type(0))) && definitely(
//						is_LT(lambda, scalar_type(1)))) {
				if (true) {
					//err_upp = rho_E_Omega(ef, eb, d, lambda);
					// computes scalar_product(math::min(t * ef, (delta - t) * eb),  math::vec_abs(l));
					scalar_type err_upp = scalar_type(0);
					duration delta_minus_t = duration(1) - lambda;
					for (size_t i = 0; i < ef.size(); ++i) {
						scalar_type t_ef = lambda * ef[i];
						scalar_type one_t_eb = delta_minus_t * eb[i];
						err_upp += std::min(lambda * ef[i], delta_minus_t * eb[i]) * d_abs[i];
					}
					//std::cout << "err_upp:" << err_upp << std::endl;
					err_low = -err_upp;

					tdi.t = t + lambda * delta;
					// we don't assign tdi.delta because we don't know the remaining switching lambdas yet
					tdi.d = direction();
					scalar_type upper_bound = (scalar_type(1) - lambda)
							* rho_Omega_beg + lambda * rho_Omega_end + err_upp;
					scalar_type lower_bound = std::max(
							l_absc_forw + lambda * l_slope_forw,
							l_absc_backw + lambda * l_slope_backw) + err_low;

//					std::cout << "inserting lambda=" << lambda << " rb "
//							<< rho_Omega_beg << " re " << rho_Omega_end
//							<< " sl_f " << l_slope_forw << " sl_b "
//							<< l_slope_backw << " errlow " << err_low
//							<< " errupp" << err_upp << " low " << lower_bound
//							<< " upp " << upper_bound << std::endl;

					tdi.itv = scalar_interval(lower_bound, upper_bound);
					seq.push_back(tdi);
				}
			}
		}
		if (lambda_sw < scalar_type(1)) {
			//		std::cout << "inserting lambda_sw="<< lambda_sw << " after " << lambda << std::endl;

			scalar_type err_upp_sw = rho_E_Omega(ef, eb, d, lambda_sw);
			scalar_type err_low_sw = -err_upp_sw;

			tdi.t = t + lambda_sw * delta;
			// we don't assign tdi.delta because we don't know the remaining switching lambdas yet
			tdi.d = direction();
			scalar_type upper_bound = (scalar_type(1) - lambda_sw)
					* rho_Omega_beg + lambda_sw * rho_Omega_end + err_upp_sw;
			scalar_type lower_bound = std::max(
					l_absc_forw + lambda_sw * l_slope_forw,
					l_absc_backw + lambda_sw * l_slope_backw) + err_low_sw;
			tdi.itv = scalar_interval(lower_bound, upper_bound);
			seq.push_back(tdi);
		}
		// sorting in time
		//std::sort(before_seq,seq.end());
		seq.sort();

		/** @todo one could remove the colinear lambda breakpoints */

		// check whether the error is within the allowed bounds:
		// err(tau) <= err_begin + (tau-t)/delta*(err_end-err_begin)
		// while we're looping, also assign delta
		timed_directional_interval_sequence::iterator it = seq.begin();
		while (!err_too_big && it != seq.end()) {
			// compute the deltas
			const duration& previous_t(it->t);
			duration& previous_delta(it->delta);
			// compute the error
			scalar_type err = it->itv.upper() - it->itv.lower();
			error_type local_err_bounds = err_begin + (it->t - t) / delta * (err_end
					- err_begin);

			err_too_big = err_too_big || definitely(
					!local_err_bounds.is_satisfied(it->itv));

//			std::cout << "local error bounds:" << local_err_bounds << " measured " << rel_abs_error<scalar_type>::measure_error(it->itv) << " itv:"<<it->itv.lower() <<","<< it->itv.upper() << " too big:" << std::boolalpha << err_too_big << std::endl;
			++it;
			if (it != seq.end()) {
				previous_delta = it->t - previous_t;
			} else {
				previous_delta = t + delta - previous_t;
			}
		}
	}

	// if error too big, split
	if (err_too_big || skip_to_smaller_delta) {
		// split the interval in two parts, first and second
		// choose the time point at which to split: t+delta_first
		duration delta_first = delta / duration(2.0); // this can be any value strictly between 0 and delta
		// the allowable error bound at the split
		error_type err_split = err_begin + delta_first / delta * (err_end
				- err_begin);

		// compute a sequence for the first part
		timed_directional_interval_sequence first_seq;
		vector_type vec_Omega_middle;
		first_seq = omega_X0_delta(d, t, delta_first, err_begin, err_split,
				vec_Omega_beg, vec_Omega_middle);

		// retrieve the support vector for the middle
		// direction at the beginning of the second interval
		direction d_second = first_seq.rbegin()->d;

		// compute a sequence for the second part
		timed_directional_interval_sequence second_seq;
		second_seq = omega_X0_delta(d_second, t + delta_first,
				delta - delta_first, err_split, err_end, vec_Omega_middle,
				vec_Omega_end);

		// fix the delta of the middle element since now we now the following
		if (second_seq.begin() != second_seq.end()) {
			first_seq.rbegin()->delta = second_seq.begin()->t - first_seq.rbegin()->t;
		}

		// return the concatenated sequences
		// GF: deactivate this since it's fairly useless; using the option for something else
//		if (my_simplify_concave) {
//			first_seq.concatenate_with_temp_and_simplify(second_seq,err_end);
//		} else {
			first_seq.concatenate_with_temp(second_seq);
//		}
		return first_seq;
	} else {
		// insert the last element
		tdi.t = t + delta;
		tdi.d = d_end;
		tdi.delta = duration(0);
		tdi.itv = scalar_interval(rho_Omega_end);
		seq.push_back(tdi);

		return seq;
	}
}


inline
void spacetime_plif::compute_rho_Omega_err(scalar_type& err_low,
		scalar_type& err_upp, const vector_type& ef, const vector_type& eb,
		const vector_type& d, const duration& lambda,
		const scalar_type& err_quad_psi_upp, const scalar_type& err_lin_psi_low) {
	scalar_type err_omega_upp = rho_E_Omega(ef, eb, d, lambda);
	err_upp = err_omega_upp + err_quad_psi_upp * lambda * lambda;
	err_low = err_lin_psi_low * lambda;
}

inline
spacetime_plif::scalar_type spacetime_plif::rho_X0(const vector_type& d) {
	if (my_aff.X0) {
		return rho(d, *my_aff.X0);
	} else
		return scalar_type(0);
}

inline
spacetime_plif::vector_type spacetime_plif::rho_vec_X0(const vector_type& d) {
	if (my_aff.X0) {
		return rho_vec(d, *my_aff.X0);
	} else
		return spacetime_plif::vector_type(domain());
}

inline
const positional_vdomain& spacetime_plif::domain() {
	return my_aff.dyn.get_A().domain();
}

inline
spacetime_plif::scalar_type spacetime_plif::rho_Ub(const vector_type& d) {
	if (has_Ub()) {
		return rho(d, *my_aff.X0);
	} else
		return scalar_type(0);
}

inline
spacetime_plif::scalar_type spacetime_plif::rho_E_Omega(const vector_type& ef,
		const vector_type& eb, const vector_type& d, const duration& lambda) {
	// computes scalar_product(math::min(t * ef, (delta - t) * eb),  math::vec_abs(l));
	scalar_type x = scalar_type(0);
	duration delta_minus_t = duration(1) - lambda;
	for (size_t i = 0; i < ef.size(); ++i) {
		scalar_type t_ef = abs(lambda * ef[i]);
		scalar_type one_t_eb = abs(delta_minus_t * eb[i]);
		x += std::min(t_ef, one_t_eb) * abs(d[i]);
	}
	return x;
}

}

#endif /* SPACETIME_PLIF_OMEGA_HPP_ */
