/*
 * spacetime_plif_helpers.hpp
 *
 *  Created on: Oct 20, 2012
 *      Author: notroot
 */

#ifndef SPACETIME_PLIF_HELPERS_HPP_
#define SPACETIME_PLIF_HELPERS_HPP_

#include "spacetime_plif.h"

#include "utility/basic_warning.h"

namespace spacetime {

inline
size_t spacetime_plif::size() const {
	return my_evolutions.size();
}

inline
spacetime_plif::scalar_type spacetime_plif::get_norm(const direction& d) {
	return d.infinity_norm();
}

inline
spacetime_plif::direction spacetime_plif::get_normed_direction(
		const direction& d) {
	scalar_type d_norm = get_norm(d);
	if (is_MEQ(d_norm, scalar_type(0))) {
		throw basic_exception(
				"trying to norm a direction vector of length zero");
	}
	return d / d_norm;
}

inline
spacetime_plif::scalar_type spacetime_plif::rho(const vector_type& l,
		const support_function_provider& S) {
	scalar_type v;
	vector_type sv;
	bool is_empty, is_bounded;
	S.compute_support(l, v, sv, is_empty, is_bounded);
	if (!is_bounded) {
		std::stringstream s;
		logger::copyfmt_to(s);
		s << "Unbounded in direction: " << l << std::endl;
		s << "Unbounded set: " << std::endl << S;
		throw basic_exception(
				"Support function evaluation requested for an unbounded set:\n"
						+ s.str());
	}
	if (is_empty)
		throw std::runtime_error(
				"Support function evaluation requested for an empty set");
	return v;
}

inline
spacetime_plif::vector_type spacetime_plif::rho_vec(const vector_type& l,
		const support_function_provider& S) {
	scalar_type v;
	vector_type sv;
	bool is_empty, is_bounded;
	if (!S.computes_support_vector())
		throw std::runtime_error(
				"Support vector requested a support_function_provider that doesn't provide a support vector");
	S.compute_support(l, v, sv, is_empty, is_bounded);
	if (!is_bounded) {
		std::stringstream s;
		s << "Unbounded in direction: " << l << std::endl;
		s << "Unbounded set: " << std::endl << S;
		throw basic_exception(
				"Support function evaluation requested for an unbounded set:\n"
						+ s.str());
	}
	if (is_empty)
		throw std::runtime_error(
				"Support function evaluation requested for an empty set");
	return sv;
}

inline
spacetime_plif::plif_type spacetime_plif::linear_interpolation(
		const timed_directional_interval_sequence& seq) {
	LOGGERSW(DEBUG7, "linear_interpolation", "linear interpolation");

//	plf_type f = plf_type::create_closed(); // lower
//	plf_type g = plf_type::create_closed(); // upper

	using plif::breakpoint;
	std::vector<breakpoint> f_pts, g_pts;
	// reserve space for at least 2x as many points as in f
	f_pts.reserve(2 * seq.size());
	g_pts.reserve(2 * seq.size());

	if (seq.begin() != seq.end()) {
		timed_directional_interval_sequence::const_iterator curr_it =
				seq.begin();
		timed_directional_interval_sequence::const_iterator next_it =
				++seq.begin();

		for (; next_it != seq.end();) {
			// convert step sequence to left limit
			// insert with one step retardation
			if (is_MEQ(curr_it->t, next_it->t)) {
//				f.insert_point(curr_it->t, curr_it->itv.lower(),
//						next_it->itv.lower(), next_it->itv.lower());
//				g.insert_point(curr_it->t, curr_it->itv.upper(),
//						next_it->itv.upper(), next_it->itv.upper());
				f_pts.push_back(
						breakpoint(curr_it->t, curr_it->itv.lower(),
								next_it->itv.lower(), next_it->itv.lower()));
				g_pts.push_back(
						breakpoint(curr_it->t, curr_it->itv.upper(),
								next_it->itv.upper(), next_it->itv.upper()));
				// skip one
				++next_it;
			} else {
//				f.insert_continuous_point(curr_it->t, curr_it->itv.lower());
//				g.insert_continuous_point(curr_it->t, curr_it->itv.upper());
				f_pts.push_back(breakpoint(curr_it->t, curr_it->itv.lower()));
				g_pts.push_back(breakpoint(curr_it->t, curr_it->itv.upper()));
			}
			curr_it = next_it;
			if (next_it != seq.end())
				++next_it;
		}
		// insert last point
		if (curr_it != seq.end()) {
//			f.insert_continuous_point(curr_it->t, curr_it->itv.lower());
//			g.insert_continuous_point(curr_it->t, curr_it->itv.upper());
			f_pts.push_back(breakpoint(curr_it->t, curr_it->itv.lower()));
			g_pts.push_back(breakpoint(curr_it->t, curr_it->itv.upper()));
		}
	}

	plf_type f = plf_type::create_closed_from_temp(f_pts);
	plf_type g = plf_type::create_closed_from_temp(g_pts);

	return plif_type(f, g);
}

inline
spacetime_plif::plif_type spacetime_plif::left_step_interpolation(
		const timed_directional_interval_sequence& seq) {
	plf_type f = plf_type::create_closed(); // lower
	plf_type g = plf_type::create_closed(); // upper

	timed_directional_interval_sequence::const_iterator it = seq.begin();
	f.insert_continuous_point(it->t, it->itv.lower());
	g.insert_continuous_point(it->t, it->itv.upper());
	scalar_type f_last = it->itv.lower();
	scalar_type g_last = it->itv.upper();
	++it;
	for (; it != seq.end(); ++it) {
		f.insert_point(it->t, f_last, it->itv.lower(), it->itv.lower());
		g.insert_point(it->t, g_last, it->itv.upper(), it->itv.upper());
		f_last = it->itv.lower();
		g_last = it->itv.upper();
	}

	return plif_type(f, g);
}

inline
math::affine_map<spacetime_plif::scalar_type> spacetime_plif::to_map(
		const matrix_type& M) {
	return math::affine_map<scalar_type>(M);
}

inline
lin_constraint_system<spacetime_plif::scalar_type> spacetime_plif::plf_to_constraints(
		const direction& d, const plf_type& f) const {
	/** For now, we don't handle unbounded plfs. */
	if (f.is_left_unbounded() || f.is_right_unbounded()) {
		throw basic_exception("can not compute constraints for unbounded PLF");
	}

	lin_constraint_system<scalar_type> lin_cons;

	// walk through each of the breakpoints of f,
	// compute the slope form
	//    f(t) = a*t + b
	// and define the linear constraint
	//    d^T x <= a*t + b,
	// which is represented in canonical constraint form by
	//    (d -a)^T (x t) - b <= 0.

	bool domain_is_singular = is_MEQ(f.get_domain().lower(),f.get_domain().upper());

	// get the new domain
	positional_vdomain dom_t = d.domain();
	dom_t.add_variable(my_time_variable);
	positional_vdomain::size_type t_pos = dom_t.pos(my_time_variable);

	// get d' = [d 0]
	direction d_t = d;
	d_t.reorder(dom_t);

	typedef plf_type::list_of_points_type list_of_points;
	const list_of_points& lop = f.get_list();

	/** treat the special case where |lop|=2 and both points are identical
	 * this should be a bug, but as a quick fix let's just deal with it.
	 * @todo find the source of this */
//	if (lop.size() == 2 && domain_is_singular && is_MEQ(lop[0].get_x(), lop[1].get_x())) {
//		if (is_MEQ(lop[0].get_y(), lop[1].get_y())) {
//			scalar_type b = lop[0].get_y();
//			lin_constraint<scalar_type> con(d_t, -b, LE);
//			lin_cons.insert(con);
//			return lin_cons;
//		} else {
//			throw basic_exception(
//					"contradicting breakpoints on singleton domain");
//		}
//	}

	unsigned int warning_superimp_count = 0;

	if (lop.size() >= 2 && !domain_is_singular) {
		for (list_of_points::const_iterator it = lop.begin();
				(it != lop.end()) && ((it + 1) != lop.end()); ++it) {
			scalar_type a, b;
			// construct slope and offset
			// since it's a concave function, there should be no discontinuities
			// @todo (to be checked)

			const scalar_type& x1 = it->get_x();
			const scalar_type& y1 = it->get_y();
			const scalar_type& x2 = (it + 1)->get_x();
			const scalar_type& y2 = (it + 1)->get_y();
			// if for some reason the x values are the same don't do anything
			if (!is_MEQ(x1, x2)) {
				a = (y2 - y1) / (x2 - x1);
				b = (y1 - a * x1);
				// assign -a as coefficient of t in d_t
				d_t[t_pos] = -a;
				// add the constraint
				lin_constraint<scalar_type> con(d_t, -b, LE);
				lin_cons.insert(con);
			} else {
				++warning_superimp_count;
			}
		}
	} else if (domain_is_singular) {
			const scalar_type& x1 = lop.front().get_x();
			const scalar_type& y1 = lop.front().get_y();
			scalar_type b;
			b = y1;
			// assign -a as coefficient of t in d_t
			d_t[t_pos] = scalar_type(0);
			lin_constraint<scalar_type> con(d_t, -b, LE);
			lin_cons.insert(con);
	} else {
		throw basic_exception(
				"can not compute linear constraints for PLF with non-singleton domain and only one breakpoint");
	}

	if (warning_superimp_count > 0) {
		using plif::operator<<;
//		std::cout << lop;
		basic_warning(__FUNCTION__,
				"piecewise linear function with "+to_string(lop.size())+" points has "+to_string(warning_superimp_count)+" occurences of points overlapping (could be numerical)",
				basic_warning::MISC);
	}

	return lin_cons;
}

inline
spacetime_plif::plf_type spacetime_plif::augment_by_error(const plf_type& f,
		const error_type& err) {
	plf_type f_shifted;

	if (err.rel() > 1e-12) {
		/** Relative error <= r <=> g <= max(f/(1+r),f*(1+r)) */
		scalar_type a = err.rel() + scalar_type(1);
		plf_type f_plus_rel_err = pointwise_maximum(a * f,
				(scalar_type(1) / a) * f);
		f_shifted = pointwise_maximum(f_plus_rel_err, f + err.abs());
	} else {
		f_shifted = f + err.abs();
	}

	return f_shifted;
}

}

#endif /* SPACETIME_PLIF_HELPERS_HPP_ */
