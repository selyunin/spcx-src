/*
 * spacetime_plif_psi.hpp
 *
 *  Created on: Oct 20, 2012
 *      Author: notroot
 */

#ifndef SPACETIME_PLIF_PSI_HPP_
#define SPACETIME_PLIF_PSI_HPP_

#include "spacetime_plif.h"

namespace spacetime {

/*
 spacetime_plif::timed_directional_interval_sequence spacetime_plif::psi_err(
 const direction& d, const time_interval& I,
 const error_type& err_begin, const error_type& err_end) {
 duration t_start = I.lower();
 duration delta = I.upper() - I.lower();
 error_type err_rate = (err_end - err_begin) / delta;
 timed_directional_interval_sequence err_seq;
 if (my_aff.U) {
 // start the error at zero and leave the slack err_begin for the error resulting from Ub
 err_seq = psi_delta(d, t_start, delta, error_type(0),
 err_end - err_begin, false);
 } else {
 // No inputs, so error is zero and the time step
 // goes over the entire interval.
 timed_directional_interval tdi;
 tdi.t = t_start;
 tdi.delta = delta;
 tdi.d = d;
 tdi.itv = scalar_interval(scalar_type(0), scalar_type(0));
 err_seq = timed_directional_interval_sequence(tdi);
 }

 //		std::cout << "err seq:" << std::endl;
 //		err_seq.print();
 //		std::cout << std::endl;

 // construct bounds on psi using the error sequence
 timed_directional_interval tdi;
 duration t = duration(0);
 scalar_interval psi_itv = scalar_interval(scalar_type(0), scalar_type(0));

 tdi.t = t;
 tdi.delta = err_seq.begin()->delta;
 tdi.d = d;
 tdi.itv = psi_itv;
 timed_directional_interval_sequence psi_seq;

 // compute the integral over the const input b
 direction b_mapped = my_aff.dyn.get_b();
 scalar_interval b_integral(0);

 timed_directional_interval_sequence::const_iterator it = err_seq.begin();
 error_type err_upp(0), err_low(0);
 for (; it != err_seq.end();) {
 // keep all info from the error sequence except the interval
 tdi.t = it->t;
 tdi.d = it->d;
 // start a new interval with nonbloated psi as left bound
 tdi.delta = it->delta;
 tdi.itv = psi_itv + b_integral;
 psi_seq.push_back(tdi);

 // now add sequence points to cover the interior of the interval [t,t+delta]
 if (true) {
 // we can incur the error err_begin plus the slack accumulated so far
 error_type err_Ub_bound = err_begin + (err_rate * it->delta
 - width(psi_itv));

 // we can't use the error from err_seq since we need Ub here.
 scalar_type ub_supp(0);
 ub_supp = rho(tdi.d, *my_Ub);
 // get the bounds on the interpolation
 timed_directional_interval_sequence err_seq_Ub;
 // start the error at zero and leave the slack err_begin for the error resulting from Ub
 err_seq_Ub = psi_delta(it->d, it->t, it->delta, error_type(0),
 err_Ub_bound, true);

 scalar_interval seq_itv(psi_itv + b_integral); // psi_itv is the solution at the discrete time points, seq_itv is the envelope of the continuous function
 for (timed_directional_interval_sequence::const_iterator jt =
 err_seq_Ub.begin(); jt != err_seq_Ub.end(); ++jt) {
 // push the bloated psi as right bound
 tdi.t = jt->t;
 tdi.delta = jt->delta;
 tdi.itv = seq_itv;
 if (jt != err_seq_Ub.begin()) {
 psi_seq.push_back(tdi);
 }
 seq_itv += scalar_interval(jt->delta * ub_supp) + jt->itv;
 }
 // push the last piece as right bound
 tdi.t = tdi.t + tdi.delta;
 tdi.delta = scalar_type(0);
 tdi.itv = seq_itv;
 psi_seq.push_back(tdi);
 }

 // integral of b
 b_mapped = phi1(it->delta) * my_aff.dyn.get_b();
 scalar_type b_supp = scalar_product(it->d, b_mapped);
 std::cout << "adding " << my_aff.dyn.get_b() << " mapped to " << b_mapped << " in direction " << it->d << " gives " << b_supp << " integral " << b_integral.lower() << std::endl;
 b_integral += scalar_interval(b_supp);

 // add the bounds on the integral of U
 scalar_type u_supp(0);
 if (my_aff.U)
 u_supp = rho(it->d, *my_aff.U);
 psi_itv += scalar_interval(it->delta * u_supp) + it->itv;

 ++it;
 }
 // add the last piece
 tdi.t = tdi.t + tdi.delta;
 tdi.d = exp_AdeltaT(tdi.delta) * tdi.d;
 tdi.delta = scalar_type(0);
 tdi.itv = psi_itv;
 psi_seq.push_back(tdi);

 //		std::cout << "psi seq:" << std::endl;
 //		psi_seq.print();
 //		std::cout << std::endl;

 return psi_seq;
 } */

inline
spacetime_plif::timed_directional_interval_sequence spacetime_plif::omega_Ub_err(
		const direction& d, const duration& t, const duration& delta,
		const error_type& err_begin, const error_type& err_end) {

	// get the sequence of error terms for Ub (use_Ub=true) with non-rate error (use_rate = false).
	timed_directional_interval_sequence seq = psi_delta(d, t, delta, err_begin,
			err_end, true, false);

	// go through the sequence and add the support of Ub
	// The resulting sequence refers to the integral over U at the end of the time interval
	// -- the value at the beginning being zero.
	if (seq.begin() != seq.end()) {
		for (timed_directional_interval_sequence::iterator it = seq.begin(); it
				!= seq.end(); ++it) {
			scalar_type ub_supp = rho(it->d, *my_Ub);
			it->itv += scalar_interval(it->delta * ub_supp);
		}
	}

	return seq;
}

 inline
spacetime_plif::timed_directional_interval_sequence spacetime_plif::psi_err(
		const direction& d, const duration& t, const duration& delta,
		const error_type& err_begin, const error_type& err_end) {
	LOGGERSW(DEBUG6,"psi_err","computing evolution of nonautonomous dynamics");

	timed_directional_interval_sequence seq;

	error_type err_omega_Ub_beg(0,0), err_omega_Ub_end(0,0);
	error_type err_psi_beg(0,0), err_psi_end(0,0);
	scalar_type err_omega_percentage(0);

	err_omega_percentage = scalar_type(0.5); // division between Omega and Psi, any value strictly between 0 and 1
	if (!my_aff.dyn.get_U()) {
		err_omega_percentage = scalar_type(1);
	}
	err_omega_Ub_beg = err_begin; // in the beginning Omega takes everything
	err_omega_Ub_end = err_end * err_omega_percentage;
	err_psi_end = err_end - err_omega_Ub_end;

	//	std::cout << "err_omega_Ub_beg: " << err_omega_Ub_beg << " err_omega_Ub_end: " << err_omega_Ub_end << " err_psi_end: " << err_psi_end << std::endl;

	timed_directional_interval_sequence omega_Ub_seq = omega_Ub_err(d,
			my_time_domain.lower(), width(my_time_domain), err_omega_Ub_beg,
			err_omega_Ub_end);

//		std::cout << "Omega_Ub_err sequence from err " << err_begin << " to "
//				<< err_end << ":" << std::endl;
//		omega_Ub_seq.print(std::cout);

	// compute psi where required in the Omega_Ub sequence
	// for each value p in Ub_err, add the points:
	// (p->t,psi(p->t))
	// (p->t+delta,psi(p->t)+p->itv)
	scalar_interval psi_itv(0.0), psi_incr(0.0);
	timed_directional_interval tdi;
	duration last_delta(0);
	for (timed_directional_interval_sequence::const_iterator it =
			omega_Ub_seq.begin(); it != omega_Ub_seq.end();) {
		// first point
		tdi = *it;
		tdi.itv = psi_itv;
		seq.push_back(tdi);
		// second point (right limit)
		tdi.t += it->delta;
		tdi.delta = duration(0);
		tdi.d = direction(); // it's just a right limit point
		tdi.itv = psi_itv + it->itv;
		seq.push_back(tdi);

		direction last_d = it->d;
		last_delta = it->delta;

		// allowed error is a linear interpolation between
		// err_omega_Ub_beg and err_omega_Ub_end
		error_type allowed_err = (it->delta) / delta * err_psi_end;
		++it;

		if (it != omega_Ub_seq.end()) {

			// compute increment of psi for the next round, using the next d
			psi_incr = psi_err_value(last_d, duration(0), last_delta,
					error_type(0,0), allowed_err);

			// update psi
			psi_itv += psi_incr;
		}
	}

//	std::cout << "psi_err sequence: " << std::endl;
//	seq.print(std::cout);

	// GF: deactivate this since it's fairly useless; using the option for something else
//	if (my_simplify_concave)
//		seq.simplify_concave(seq.begin(),seq.end());

	return seq;
}

 inline
spacetime_plif::scalar_interval spacetime_plif::psi_err_value(
		const direction& d, const duration& t, const duration& delta,
		const error_type& err_begin, const error_type& err_end) {
	// returns the value of psi at exactly time t+delta

	// if delta is zero, psi is also zero
	if (is_MEQ(delta, duration(0))) {
		return scalar_interval(scalar_type(0));
	}

	error_type err_rate = (err_end - err_begin) / delta;

	timed_directional_interval_sequence err_seq;
	if (my_aff.dyn.get_U()) {
		// start the error at zero and leave the slack err_begin for the error resulting from Ub
		err_seq = psi_delta(d, t, delta, error_type(0,0), err_end - err_begin,
				false, true);
	} else {
		// No inputs, so error is zero and the time step
		// goes over the entire interval.
		timed_directional_interval tdi;
		tdi.t = t;
		tdi.delta = delta;
		tdi.d = d;
		tdi.itv = scalar_interval(scalar_type(0), scalar_type(0));
		err_seq = timed_directional_interval_sequence(tdi);
	}

	//		std::cout << "err seq:" << std::endl;
	//		err_seq.print();
	//		std::cout << std::endl;

	// construct bounds on psi using the error sequence
	timed_directional_interval tdi;
	scalar_interval psi_itv = scalar_interval(scalar_type(0), scalar_type(0));

	tdi.t = t;
	tdi.delta = err_seq.begin()->delta;
	tdi.d = d;
	tdi.itv = psi_itv;
	timed_directional_interval_sequence psi_seq;

	// compute the integral over the const input b
	direction b_mapped = my_aff.dyn.get_b();
	scalar_interval b_integral(0);

	timed_directional_interval_sequence::const_iterator it = err_seq.begin();
	scalar_type err_upp(0), err_low(0);
	for (; it != err_seq.end();) {
		/* Don't need tdi here
		 // keep all info from the error sequence except the interval
		 tdi.t = it->t;
		 tdi.d = it->d;
		 // start a new interval with nonbloated psi as left bound
		 tdi.delta = it->delta;
		 tdi.itv = psi_itv + b_integral;
		 psi_seq.push_back(tdi);
		 */

		// integral of b
		b_mapped = phi1(it->delta) * my_aff.dyn.get_b();
		scalar_type b_supp = scalar_product(it->d, b_mapped);
		//		std::cout << "adding " << my_aff.dyn.get_b() << " mapped to "
		//				<< b_mapped << " in direction " << it->d << " gives " << b_supp
		//				<< " integral " << b_integral.lower() << std::endl;
		b_integral += scalar_interval(b_supp);

		// add the bounds on the integral of U
		scalar_type u_supp(0);
		if (my_aff.dyn.get_U())
			u_supp = rho(it->d, *my_aff.dyn.get_U());
		psi_itv += scalar_interval(it->delta * u_supp) + it->itv;

		++it;
	}
	// add the last piece
	/* no need for the sequence, we just want psi_tv
	 tdi.t = tdi.t + tdi.delta;
	 tdi.d = exp_AdeltaT(tdi.delta) * tdi.d;
	 tdi.delta = scalar_type(0);
	 tdi.itv = psi_itv + b_integral;
	 psi_seq.push_back(tdi);
	 */

	//	std::cout << "psi seq for err value:" << std::endl;
	//	psi_seq.print();
	//	std::cout << std::endl;

	//return psi_seq;

	// return the values at the end of the interval
	return psi_itv + b_integral;
}

 inline
spacetime_plif::timed_directional_interval_sequence spacetime_plif::psi_delta(
		const direction& d, const duration& t, const duration& delta,
		const error_type& err_begin, const error_type& err_end, bool use_Ub,
		bool use_rate) {

	scalar_type err_upp;
	scalar_type err_low;

	bool skip_to_smaller_delta = (min_infeasible_delta > scalar_type(0) && delta >= min_infeasible_delta/scalar_type(16));

	if (!skip_to_smaller_delta) {
		// try to obtain the parameters for *this; if there are numerical issues
		// skip to next recursion with smaller delta
		try {
			// compute the upper error bound
			if (use_Ub) {
				err_upp = rho(d, E_Ub(delta));
			} else {
				err_upp = rho(d, E_U(delta));
			}

			// compute the lower error bound
			// GF: added err_upp since otherwise it's not conservative
			if (use_Ub) {
				err_low = -err_upp - rho(d, *negAPhi2Ub(delta));
			} else {
				err_low = -err_upp - rho(d, *negAPhi2U(delta));
			}
		} catch (math::invalid_number_exception& e) {
//		std::cout << std::endl
//				<< "decreasing delta in psi because of numerical issues"
//				<< std::endl;
			LOGGER(DEBUG6, __FUNCTION__,
					"decreasing current time step "+to_string(delta)+" because of numerical issues");
			skip_to_smaller_delta = true;
		}
	}
	// derive the error bound depending on the mode
	error_type allowed_error;
	bool error_too_big = false;
	if (!skip_to_smaller_delta) {
		if (use_rate) {
			// check whether the error conforms to the error rate
			allowed_error = err_end - err_begin;
			error_too_big = !maybe(allowed_error.is_satisfied(err_low,err_upp));
		} else {
			// check whether the error lies below the lowest
			// allowed_error = std::min(err_end, err_begin);
			error_too_big = !(maybe(err_begin.is_satisfied(err_low,err_upp)) && maybe(err_end.is_satisfied(err_low,err_upp)));
		}
	}
	if (skip_to_smaller_delta || error_too_big) {
		// split the interval in two parts, first and second
		// choose the time point at which to split: t+delta_first
		duration delta_first = delta / duration(2.0); // this can be any value strictly between 0 and delta

		// compute a sequence for the first part
		timed_directional_interval_sequence first_seq;
		first_seq = psi_delta(d, t, delta_first, err_begin,
				err_begin + (err_end - err_begin) * delta_first / delta,
				use_Ub, use_rate);

		/** @todo Use as error start the error accumulated in first_seq */

		// compute a sequence for the second part
		timed_directional_interval_sequence second_seq;
		// direction at the beginning of the second interval

		// The following may throw because of numerics
		//		direction d_second = exp_AdeltaT(delta_first) * d;
		// compute it from the last element in first_seq so that delta is not
		// too big (numerical issues)
		direction d_middle = first_seq.rbegin()->d;
		duration t_middle = first_seq.rbegin()->t;
		duration t_second = t + delta_first;
		duration delta_second = delta - delta_first;
		duration delta_middle = t_second-t_middle;
		direction d_second = d_middle;
		while (definitely(is_LT(t_middle, t_second))) {
			try {
				d_second = exp_AdeltaT(delta_middle) * d_second;
				t_middle += delta_middle;
			} catch (math::invalid_number_exception& e) {
				delta_middle = delta_middle / duration(2.0);
				LOGGER(DEBUG6,__FUNCTION__,"decreasing intermediate time step to "+to_string(delta_middle)+" because of numerical issues");
			}
		}
		if (definitely(is_GT(t_middle, t_second))) {
			throw basic_exception("something went wrong with computing time steps in "+std::string(__FUNCTION__));
		}
		second_seq = psi_delta(d_second, t_second, delta_second,
				err_begin + (err_end - err_begin) * delta_first / delta,
				err_end, use_Ub, use_rate);

		// return the concatenated sequences
		first_seq.concatenate_with_temp(second_seq);

		return first_seq;
	} else {
		// no need to refine the sequence, return single point
		timed_directional_interval tdi;
		tdi.t = t;
		tdi.delta = delta;
		tdi.d = d;
		tdi.itv = scalar_interval(err_low, err_upp);
		timed_directional_interval_sequence seq(tdi);
		return seq;
	}
}

}

#endif /* SPACETIME_PLIF_PSI_HPP_ */
