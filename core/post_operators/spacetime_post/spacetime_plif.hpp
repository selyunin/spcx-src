/*
 * spacetime_plif.hpp
 *
 *  Created on: Oct 20, 2012
 *      Author: notroot
 */

#ifndef SPACETIME_PLIF_HPP_
#define SPACETIME_PLIF_HPP_

#include "spacetime_plif.h"
#include "core/continuous/support_function/sf_derived/sf_sum.h"

#include "math/vdom/lin_constraint_system.h"
#include "extern/plif/piecewise_linear_function_operators.h" // just for temp debug
#include "extern/plif/piecewise_linear_interval_function_operators.h"

namespace spacetime {

inline spacetime_plif::spacetime_plif(
		affine_support_flowpipe_description<scalar_type> prob,
		time_interval intv) :
		my_aff(prob), my_time_domain(intv) {
	// consistency check
	if (my_aff.dyn.domain() != my_aff.dyn.codomain()) {
		throw std::runtime_error("inconsistent dynamics: nonsquare matrix");
	}

	// initialize the delta_parameters pointer
	my_delta_params_it = my_delta_params_cache.end();
	min_infeasible_delta = scalar_type(-1);
	min_delta = scalar_type(-1);

	// compute parameters that are independent of delta
	const matrix_type& A = my_aff.dyn.get_A();
	if (my_aff.X0) {
		affine_map<scalar_type> A2_map = to_map(A * A);
		my_S2f = finite_symmetric_bounding_box(*my_aff.X0, A2_map);
	}

	// compute Ub
	my_Ub = my_aff.dyn.get_U(); // start with U+b=U in case there isn't any nonzero b
	support_function_provider::ptr b_box;
	// if there is a nonzero b, create the corresponding sets
	bool b_zero = is_MEQ(my_aff.dyn.get_b().infinity_norm(), scalar_type(0));
	//std::cout << std::endl << my_aff.dyn.get_b().get_vector() << std::endl;
	if (my_aff.dyn.get_b().size() > 0 && !b_zero) {
		const positional_vdomain& dom = A.domain();
		vector_type zeros(dom);

		// create a support function provider representing b
		b_box = support_function_provider::ptr(
				new finite_hyperbox<scalar_type>(my_aff.dyn.get_b(), zeros));

		// if there is U and nonzero b, put the two together
		if (my_aff.dyn.get_U()) {
			my_Ub = support_function_provider::ptr(
					new support_function::sf_sum<scalar_type>(my_aff.dyn.get_U(), b_box));
		} else {
			my_Ub = b_box;
		}
	}

	// get a time variable
	// for now: always use the same
	my_time_variable = variable("SPACETIME_TIME");

	// default cut point method
	my_cut_point_method.type = cut_point_method::FIXED_COUNT_SAME_SIZE_PIECES;
	my_cut_point_method.piece_count = 1;

	// by default, don't simplify concanve
	my_simplify_concave = false;
	my_simplify_convex = false;
}

inline const positional_vdomain& spacetime_plif::domain() const {
	return my_aff.dyn.domain();
}

inline const spacetime_plif::time_interval& spacetime_plif::get_time_domain() const {
	return my_time_domain;
}

inline const affine_support_flowpipe_description<spacetime_plif::scalar_type>& spacetime_plif::get_affine_support_flowpipe_description() const {
	return my_aff;
}

inline
void spacetime_plif::set_simplify_concave(bool b) {
	my_simplify_concave = b;
}

inline
void spacetime_plif::set_simplify_convex(bool b) {
	my_simplify_convex = b;
}

inline const spacetime_plif::annotated_plif& spacetime_plif::get_or_compute_evolution(
		direction d, const error_type& err) {
	d.reorder(domain());
	// consistency check
	if (domain().size() != d.size()) {
		throw std::runtime_error(
				"trying to compute evolution in dimension "
						+ to_string(my_aff.dyn.get_A().domain().size())
						+ " for a direction of dimension "
						+ to_string(d.size()));
	}

	// check if d is already in the store
	direction d_normed = get_normed_direction(d);
	evolution_cache::iterator it = my_evolutions.find(d_normed);
	if (it == my_evolutions.end()) {
		// direction not found, compute it
		annotated_plif new_plif = compute_evolution(d_normed, err);
		it = my_evolutions.insert_missing(d_normed, new_plif);
	} else if (!maybe(it->second.second <= err)) {
		// recompute because the error is too big
		annotated_plif new_plif = compute_evolution(d_normed, err);
		// update cache with more precise result
		it->second = new_plif;
	}
	return it->second;
}

inline const spacetime_plif::annotated_plif* spacetime_plif::get_evolution(
		direction d) {
	d.reorder(domain());

	// check if d is already in the store
	direction d_normed = get_normed_direction(d);
	evolution_cache::iterator it = my_evolutions.find(d_normed);
	if (it == my_evolutions.end()) {
		return 0;
	} else {
		return &it->second;
	}
}

inline spacetime_plif::annotated_plif spacetime_plif::compute_evolution(
		const direction& d, const error_type& err) {
	LOGGERSW(DEBUG7, "compute_evolution", "computing support evolutions");
	logger::logger_id sw_id = logger::get_last_id();
	logger::logger_id log_id;
	IFLOGGER(DEBUG7) {
		log_id =
				LOGGER(DEBUG5,"compute_evolution","computing evolution in direction "+logger::formatted_to_string(d)+" up to err "+to_string(err));
	}
	plif_type evo;

	//std::cout << "CE time domain: " << my_time_domain << ", width: " << width(my_time_domain);
	if (math::definitely(is_LT(width(my_time_domain),scalar_type(0)))) {
		throw std::runtime_error("negative duration");
	}

	/** Divide the error margin up between omega_X0, omega_Ub and psi
	 *
	 * Since the error of psi is forcibly increasing over time, we adopt
	 * a linear model:
	 * Initially, the error of psi is zero, so omega_X0 and omega_Ub can
	 * use all the margin. Omega_X0 takes up X0_percentag and Omega_Psi the rest.
	 * At the end, omega_X0 and omega_ub must share the error margin with psi.
	 * They take up Omega_percentage, and psi the rest.
	 *
	 */
	timed_directional_interval_sequence omega_X0_seq;
	scalar_type err_omega_X0_beg(0), err_omega_X0_end(0);
	scalar_type err_omega_Ub_beg(0), err_omega_Ub_end(0);
	scalar_type err_psi_beg(0), err_psi_end(0);
	scalar_type err_omega_percentage(0), err_X0_percentage(0);

	// the relative error is reserved for Omega_X0, the absolute error is divided
	// equally between Omega_X0 and all that's Psi

	err_omega_percentage = scalar_type(0.5); // division between Omega and Psi, any value strictly between 0 and 1
	if (my_aff.X0) {
		if (has_Ub()) {
			err_X0_percentage = scalar_type(0.5); // division between Omega_X0 and Omega_Ub, any value strictly between 0 and 1
		} else { // no U, so X can take all the error margin
			err_omega_percentage = scalar_type(1); // Psi is zero, so Omega can use everthing
			err_X0_percentage = scalar_type(1); // Ub is zero, so X0 can use everything
		}
	} else { // no X0, we can assume there's Ub
		err_X0_percentage = scalar_type(0);
	}
	err_omega_X0_beg = err.abs() * err_X0_percentage; // in the beginning Omega takes everything
	err_omega_X0_end = err.abs() * err_omega_percentage * err_X0_percentage;
	err_omega_Ub_beg = err.abs() - err_omega_X0_beg; // in the beginning Omega takes everything
	err_omega_Ub_end = err.abs() * err_omega_percentage
			* (scalar_type(1) - err_X0_percentage);
	err_psi_end = err.abs() - err_omega_X0_end - err_omega_Ub_end;

	// compute the omega part coming from X0
	if (my_aff.X0) {
		if (!my_simplify_concave) {
			omega_X0_seq = omega_X0_err(d, my_time_domain.lower(),
					width(my_time_domain),
					error_type(err.rel(), err_omega_X0_beg),
					error_type(err.rel(), err_omega_X0_end));
		} else {
			omega_X0_seq = omega_X0_err_simplified(d, my_time_domain.lower(),
					width(my_time_domain),
					error_type(err.rel(), err_omega_X0_beg),
					error_type(err.rel(), err_omega_X0_end));
		}
		LOGGER(DEBUG6, "compute_evolution",
				"computed omega sequence of size "+to_string(omega_X0_seq.size()));
	}

	// compute the omega part coming from Ub
	timed_directional_interval_sequence psi_seq;
	if (has_Ub()) {
		psi_seq = psi_err(d, my_time_domain.lower(), width(my_time_domain),
				error_type(scalar_type(0), err_omega_Ub_beg),
				error_type(scalar_type(0), err_omega_Ub_end + err_psi_end));
		LOGGER(DEBUG6, "compute_evolution",
				"computed psi sequence of size "+to_string(psi_seq.size()));
	}

//	std::cout << std::endl << "Omega: "<< std::endl; omega_X0_seq.print(std::cout);
//	std::cout << std::endl << "Psi: "<< std::endl; psi_seq.print(std::cout);

	plif_type res_plif;
	// Convert to plifs and add
	if (my_aff.X0) {
		if (has_Ub()) {
			LOGGERSW(DEBUG6, "compute_evolution", "superposing");

			plif_type X0_plif = linear_interpolation(omega_X0_seq);
			plif_type psi_plif = linear_interpolation(psi_seq);
			res_plif = X0_plif + psi_plif;
		} else {
			res_plif = linear_interpolation(omega_X0_seq);
		}
	} else if (has_Ub()) {
		res_plif = linear_interpolation(psi_seq);
	}

	evo = plif::simplify(res_plif);

//	std::cout << "final evolution: " << std::endl;
//	res_plif.display();
	assert(plif::indistinguishable_count(evo.get_upper())==0);

	//error_type err_max(scalar_type(0),max_width(evo));
	annotated_plif h(evo, err);

//	IFLOGGER(DEBUG5) {
//		LOGGER_ATTACH(DEBUG5,"compute_evolution",", measured err "+to_string(max_width(evo)),log_id);
//	}
	LOGGER_ATTACH(DEBUG7, "compute_evolution",
			" with "+to_string(res_plif.size())+" points", sw_id);
	return h;
}

inline void spacetime_plif::assign_convex_hull() {
	// For each of the evolutions, replace the upper bound by its concave hull
	for (evolution_cache::iterator it = my_evolutions.begin();
			it != my_evolutions.end(); ++it) {
			const plf_type& upper_plf = it->second.first.get_upper();
			plf_type chull = concave_hull(upper_plf);

//			plif::plf_graph(chull,"gif","/tmp/test_greedy");
//			plif::piecewise_linear_interval_function margin(it->second.first.get_lower(),it->second.first.get_upper());
//			plif::plif_graph(margin,"X","/tmp/test_spacetime_margin_greedy","-m3 -W 0.006 -S 4 /tmp/test_greedy.txt -s -W 0 ");

			it->second.first.set_upper(chull);
			// @todo replace error by measured error
	}
}

inline
void spacetime_plif::set_cut_point_method(cut_point_method m) {
	my_cut_point_method = m;
}

inline spacetime_plif::cut_point_method spacetime_plif::get_cut_point_method() const {
	return my_cut_point_method;
}

inline const variable& spacetime_plif::get_time_variable() const {
	return my_time_variable;
}

inline std::vector<spacetime_plif::time_interval> spacetime_plif::get_maybe_satisfying_subdomains(
		const lin_constraint<scalar_type>& con, const error_type& err, bool widen_by_error) {
	using plif::operator<<;

	std::vector<spacetime_plif::time_interval> sat_intervals;
	lin_constraint<scalar_type> c = con.get_canonical();
	bool is_nonempty = true;

	if (con.is_equality()) {
		// if equality, restrict with LE first, then with GE
		c.set_sign(LE);
		sat_intervals = get_maybe_satisfying_subdomains(c, err, widen_by_error);

//		std::cout << "get_maybe_satisfying_subdomains EQ1:" << c << " sat when "<< sat_intervals << std::endl;
//		LOGGER_OS(DEBUG7, __FUNCTION__) << "LE sat domains: " << sat_intervals;

		c = lin_constraint<scalar_type>(-c.get_l(), LE);
		if (!sat_intervals.empty()) {
			std::vector<spacetime_plif::time_interval> sat_intervals2 =
					get_maybe_satisfying_subdomains(c, err, widen_by_error);
//			LOGGER_OS(DEBUG7, __FUNCTION__) << "GE sat domains: " << sat_intervals2;
			sat_intervals = plif::intersect_intervals(sat_intervals,
					sat_intervals2);
//			std::cout << "get_maybe_satisfying_subdomains EQ2:" << c << " sat when "<< sat_intervals2 << " intersected " << sat_intervals << std::endl;
		}
	} else {
		// get the evolution that corresponds to the constraint
		const direction& a = c.get_normal();
		scalar_type b = -c.get_canonic_inh_coeff(); // for constraint a^T x <= b
		// a^T x <= b is the same as -a^T x >= -b
		//evo=get_or_compute_evolution(d,err)*get_norm(d);
		spacetime_plif::annotated_plif evo_ann = get_or_compute_evolution(-a,
				err);

		plf_type upper_bound;
		if (widen_by_error) {
			upper_bound = augment_by_error(evo_ann.first.get_lower(),err);
		} else {
			upper_bound = evo_ann.first.get_upper();
		}

		sat_intervals = plif::above_threshold_subdomains(
				upper_bound, -b/get_norm(a));
//		std::cout << "get_maybe_satisfying_subdomains LE:" << a << "<=" << b << " sat when "<< sat_intervals << std::endl;
	}
	return sat_intervals;
}

inline std::vector<spacetime_plif::time_interval> spacetime_plif::get_definitely_satisfying_subdomains(
		const lin_constraint<scalar_type>& con, const error_type& err) {
	using plif::operator<<;

	std::vector<spacetime_plif::time_interval> sat_intervals;
	lin_constraint<scalar_type> c = con.get_canonical();
	bool is_nonempty = true;

	if (con.is_equality()) {
		// if equality, restrict with LE first, then with GE
		c.set_sign(LE);
		sat_intervals = get_definitely_satisfying_subdomains(c, err);

//		std::cout << "get_maybe_satisfying_subdomains EQ1:" << c << " sat when "<< sat_intervals << std::endl;

		c = lin_constraint<scalar_type>(-c.get_l(), LE);
		if (!sat_intervals.empty()) {
			std::vector<spacetime_plif::time_interval> sat_intervals2 =
					get_definitely_satisfying_subdomains(c, err);
//			std::cout << "get_maybe_satisfying_subdomains EQ2:" << c << " sat when "<< sat_intervals2 << std::endl;
			sat_intervals = plif::intersect_intervals(sat_intervals,
					sat_intervals2);
		}
	} else {
		// get the evolution that corresponds to the constraint
		const direction& a = c.get_normal();
		scalar_type b = -c.get_canonic_inh_coeff(); // for constraint a^T x <= b
		// max a^T x <= b is the same as - max a^T x >= -b
		//evo=get_or_compute_evolution(d,err)*get_norm(d);
		spacetime_plif::annotated_plif evo_ann = get_or_compute_evolution(a,
				err);
		sat_intervals = plif::above_threshold_subdomains(
				-evo_ann.first.get_upper(), -b/get_norm(a));
//		std::cout << "get_maybe_satisfying_subdomains LE:" << a << "<=" << b << " sat when "<< sat_intervals << std::endl;
	}
	return sat_intervals;
}

inline
void spacetime_plif::restrict_to_constraints(const lin_constraint_system<scalar_type>& cons) {
//	using namespace plif;
	using plif::breakpoint;
	if (!cons.empty()) {
		// use a polyhedron to compute the support in cons
		typedef continuous::constr_polyhedron<scalar_type> poly_type;
		poly_type poly;
		poly.add_constraints(cons);

		// For each of the N evolutions, restrict the plif to the max induced by the constraints
		for (evolution_cache::iterator it = my_evolutions.begin();
				it != my_evolutions.end(); ++it) {
			// check if limited by con
			const direction& d = it->first;
			plif_type& plif = it->second.first;

			// compute the support
			scalar_type max_value;
			direction support_vec;
			bool is_empty,is_bounded;
			poly.compute_support(d,max_value,support_vec,is_empty,is_bounded);

			if (is_bounded) {
				// check if con is active on *it
				scalar_type plif_max = plif::supremum(plif.get_upper());
				if (definitely(is_LT(max_value, plif_max))) {
					const plf_type& f = plif.get_upper();
					// restrict upper bound of plif
					// take the first and the last point of the plif and use the max to obtain interpolation
					plf_type::list_of_points_type pts;
					pts.push_back(
							breakpoint(f.get_domain().lower(), max_value));
					pts.push_back(
							breakpoint(f.get_domain().upper(), max_value));
					plf_type g = plf_type::create_closed_from_temp(pts);
					plif.set_upper(pointwise_minimum(f, g));
				}
			}
		}
	}
}

inline
void spacetime_plif::restrict_to_subdomain(const time_interval& intv) {
	const time_interval& old_domain = get_time_domain();

	// don't do anything if intv old_domain contains intv
	if (math::definitely(is_GT(intv.lower(), old_domain.lower()))
			|| math::definitely(is_LT(intv.upper(), old_domain.upper()))) {
		time_interval new_domain = plif::intersect(get_time_domain(), intv);
		LOGGER(DEBUG7, __FUNCTION__, "restricting domain to "+to_string(new_domain.upper()));

		using plif::operator<<;
//std::cout << "restricting " << get_time_domain() << " to " << intv << " gives " << new_domain << std::endl;

		if (!plif::empty(new_domain)) {
			// @to do: adapt the problem description so it fits (future evolutions!)

			// we need a vector of intervals (even if it's just one)
			std::vector<time_interval> intv_vec(1, new_domain);

			// For each of the N evolutions, restrict the plif to the new domain
			for (evolution_cache::iterator it = my_evolutions.begin();
					it != my_evolutions.end(); ++it) {
				std::vector<plif_type> plif_pieces = dissect(it->second.first,
						intv_vec);

				//std::cout << std::endl << "pieces: " << std::endl;
				//plif_pieces.begin()->display(); //->get_upper().display();
				//(++plif_pieces.begin())->get_upper().display();

				it->second.first = *plif_pieces.begin();
			}

			// fix the state variables
			// @todo
		} else {
			throw std::runtime_error("restricting flowpipe to empty domain, don't know what to do");
		}

		my_time_domain = new_domain;
	}
}

inline math::tribool spacetime_plif::decide_outer_contains(spacetime_plif& f,
		bool refine) const {
	LOGGERSWOC(DEBUG, "decide_outer_contains", "deciding containment of flowpipes");
	LOGGER(DEBUG6, "decide_outer_contains", "deciding containment of flowpipes");
	/** @attention This is a hack to avoid a bug in stc scenario that
	 * results from putting a state on the passed as well as the waiting list.
	 * In the merging process this state will be removed from the waiting.
	 */
	if (&f==this) return math::indeterminate();

	math::tribool contains_f = true; // start optimistic

	// compare existing flowpipe evolutions
	plif_type common_time_map;
	for (evolution_cache::const_iterator it = my_evolutions.begin();
			it != my_evolutions.end(); ++it) {

		const direction& d = it->first;
		const annotated_plif& evo_ann = it->second;
		const annotated_plif* evo_ann_f = f.get_evolution(d);
		if (!evo_ann_f) {
			LOGGER(DEBUG5, __FUNCTION__, "couldn't find evolution, trying to refine");
			if (refine) {
				evo_ann_f = &f.get_or_compute_evolution(d, evo_ann.second);
			} else {
				return math::indeterminate();
			}
		}
		plif_type time_map = compute_inclusion_map(evo_ann_f->first.get_upper(),
				evo_ann.first.get_upper());

		if (plif::empty(time_map.get_domain())) {
			LOGGER(DEBUG7, __FUNCTION__, "inclusion map empty");
			return false;
		}
		if (math::definitely(is_GT(time_map.get_domain().lower(), f.get_time_domain().lower()))
				|| math::definitely(
						is_LT(time_map.get_domain().upper(), f.get_time_domain().upper()))) {
//			using plif::operator<<;
//			std::cout << "map: " << time_map << ", domain: "<<f.get_time_domain()<< std::endl;
			LOGGER(DEBUG7, __FUNCTION__, "inclusion map doesn't cover whole time domain");
			return false;
		}

		//plif_graph(time_map, "X", "/tmp/containment");
		if (it == my_evolutions.begin())
			common_time_map = time_map;
		else
			common_time_map = intersection(common_time_map, time_map);
		//plif_graph(common_time_map, "X", "/tmp/containment");
	}

	// plif_graph(common_time_map, "X", "/tmp/containment");

	// Check first if the time domain spans everything
	if (plif::empty(common_time_map.get_domain())) {
		LOGGER(DEBUG7, __FUNCTION__, "inclusion map empty");
		return false;
	}
	if (math::definitely(
			is_GT(common_time_map.get_domain().lower(),
					f.get_time_domain().lower()))
			|| math::definitely(
					is_LT(common_time_map.get_domain().upper(),
							f.get_time_domain().upper()))) {
		LOGGER(DEBUG7, __FUNCTION__, "inclusion map doesn't cover whole time domain");
		return false;
	}

	// test containment by verifying that the map is connected
	if (!plif::empty(common_time_map.get_domain())) {
		plf_type e = common_time_map.get_upper() - common_time_map.get_lower();
		double diff = infimum(e); // if negative, lower is higher than upper
		contains_f = math::maybe(is_GE(diff, 0.0));
		LOGGER(DEBUG7, __FUNCTION__, "inclusion map computes e="+to_string(diff));
	} else {
		contains_f = false;
	}

	// @todo If contains_f is false, test whether lower bounds are contained. If yes,
	// recompute upper bounds with lower error.

	return contains_f;
}

inline void spacetime_plif::set_numeric_error(error_type e) {
	my_numeric_error = e;
	plif::numeric_less<scalar_type>::set_abs_err(e.abs());
	plif::numeric_less<scalar_type>::set_rel_err(e.rel());
}

inline
support_function_provider::ptr spacetime_plif::get_last_states() const {
	support_function_provider::ptr s;
	using namespace continuous::support_function;

	double T = get_time_domain().upper();

	if (true) {
	// the following causes numerical troubles, amongst other things
//	// get e^AT*X0 + b_mapped
//	matrix_type Phi = get_delta_params(T).Phi;
//	vector_type b_mapped = phi1(T) * my_aff.dyn.get_b();
//
//	math::affine_map<double> M(Phi, b_mapped);
//	sf_unary<double>::ptr X0_mapped;
//	if (my_aff.X0) {
//		X0_mapped = sf_unary<double>::ptr(new sf_unary<double>(my_aff.X0, M));
//		s = X0_mapped;
//	}
//	if (my_aff.U) {
		// for now, just take outer approx

		lin_constraint_system<scalar_type> lin_cons;
		for (evolution_cache::const_iterator it =
				my_evolutions.begin(); it != my_evolutions.end(); ++it) {

			const direction& d = it->first;
			const annotated_plif& evo_ann = it->second;
			const plf_type& f = evo_ann.first.get_upper();

			// get U at last point
			// to do so, take outer evolution and subtract X0
			typedef plf_type::list_of_points_type list_of_points;
			const list_of_points& lop = f.get_list();

			const scalar_type& y = lop.back().get_y();
			scalar_type b = y;
			if (s) {
				b -= rho(d, *s);
			} else {
				// don't substract b_mapped
				// (otherwise we'd just add it back in at the end anyway)
			}
			lin_constraint<scalar_type> con(d, -b, LE);
			lin_cons.insert(con);
		}
		constr_polyhedron<scalar_type>::ptr poly(
				new constr_polyhedron<scalar_type>());
		poly->add_constraints(lin_cons);
		if (s) {
			s = sf_sum<double>::ptr(new sf_sum<double>(s,poly));
		} else {
			s = poly;
		}
	}

	return s;
}

inline
void spacetime_plif::print(std::ostream& os) const {
	for (evolution_cache::const_iterator it =
			my_evolutions.begin(); it != my_evolutions.end(); ++it) {

		const direction& d = it->first;
		const annotated_plif& evo_ann = it->second;

		os << "direction: " << d << ", plif: " << evo_ann.first << ", error: " << evo_ann.second << std::endl;

	}
}

}

#endif /* SPACETIME_PLIF_HPP_ */
