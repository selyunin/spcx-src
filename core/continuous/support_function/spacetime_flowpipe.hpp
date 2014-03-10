/*
 * spacetime_flowpipe.hpp
 *
 *  Created on: Nov 3, 2012
 *      Author: notroot
 */

#ifndef SPACETIME_FLOWPIPE_HPP_
#define SPACETIME_FLOWPIPE_HPP_

#include "spacetime_flowpipe.h"

#include "core/post_operators/spacetime_post/spacetime_plif.h"

namespace continuous {

template<typename scalar_type>
spacetime_flowpipe<scalar_type>::spacetime_flowpipe(
		const spacetime_flowpipe<scalar_type>& flp) : index_to_variable_id_map_provider(flp),my_splif(flp.my_splif) {
	if (flp.my_map)
	my_map = affine_map_ptr(new affine_map(*flp.my_map));
	my_W = flp.my_W;
	my_refinement_locked = flp.my_refinement_locked;
	my_outer_constraints=flp.my_outer_constraints;
}


template<typename scalar_type>
const typename spacetime_flowpipe<scalar_type>::time_interval& spacetime_flowpipe<
		scalar_type>::get_time_domain() const {
	return my_splif.get_time_domain();
}

template<typename scalar_type>
std::vector<spacetime_flowpipe<scalar_type> > spacetime_flowpipe<scalar_type>::convexify(
		const cut_point_method& m) const {
	std::vector<spacetime_plif> splifs = my_splif.convexify(m);

	std::vector<spacetime_flowpipe<scalar_type> > convex_flowpipes;
	convex_flowpipes.reserve(splifs.size());

	// copy the state
	for (size_t i = 0; i < splifs.size(); ++i) {
		convex_flowpipes.push_back(
				spacetime_flowpipe<scalar_type>(splifs[i], my_map, my_W,
						my_refinement_locked, my_outer_constraints));
	}
	return convex_flowpipes;
}

template<typename scalar_type>
polyhedron_collection<scalar_type> spacetime_flowpipe<scalar_type>::compute_outer_polyhedra(
		const cut_point_method& m) const {
	polyhedron_collection<scalar_type> polys;
	polys = my_splif.compute_outer_polyhedra(m);
	polys.add_constraints(my_outer_constraints);
	return polys;
}

template<typename scalar_type>
typename polyhedron<scalar_type>::ptr spacetime_flowpipe<scalar_type>::compute_convex_outer_polyhedron(error_type err) const {
	cut_point_method m;
	m.type = cut_point_method::FIXED_COUNT_SAME_SIZE_PIECES;
	m.piece_count = 1;
	m.approx_error = err;
	polyhedron_collection<scalar_type> polys = compute_outer_polyhedra(m);
	assert(polys.size() == 1);

	return *polys.begin();
}

template<typename scalar_type>
void spacetime_flowpipe<scalar_type>::add_outer_constraint(
		const math::lin_constraint<scalar_type>& con, bool restrict_evolutions) {
	my_outer_constraints.insert(con);
	if (restrict_evolutions) {
		math::lin_constraint_system<scalar_type> cons;
		cons.insert(con);
		my_splif.restrict_to_constraints(cons);
	}
}

template<typename scalar_type>
void spacetime_flowpipe<scalar_type>::add_outer_constraints(
		const math::lin_constraint_system<scalar_type>& cons, bool restrict_evolutions) {
	for (typename math::lin_constraint_system<scalar_type>::const_iterator it =
			cons.begin(); it != cons.end(); ++it) {
		// don't yet restrict evolutions, we'll do it later for all constraints together
		add_outer_constraint(*it,false);
	}
	if (restrict_evolutions) {
		my_splif.restrict_to_constraints(cons);
	}
}

template<typename scalar_type>
void spacetime_flowpipe<scalar_type>::restrict_to_maybe_sat_prefix(
		const math::lin_constraint_system<scalar_type>& cons,
		const error_type& err) {
	using namespace math::numeric;

	LOGGERSW(DEBUG4, "restrict_to_maybe_sat_prefix",
			"restricting to maybe satisfying prefix");

	const time_interval& tdom(my_splif.get_time_domain());
	time_interval sat_prefix(my_splif.get_time_domain());
	// get the prefix interval on which a constraint is satisfied
	for (typename math::lin_constraint_system<scalar_type>::const_iterator it =
			cons.begin();
			it != cons.end() && !plif::empty(sat_prefix); ++it) {
		// remap to domain if necessary
		math::lin_expression<scalar_type> expr = it->get_l();
		expr.reorder(domain());
		math::lin_constraint<scalar_type> con(expr, it->get_sign());
		IFLOGGER(DEBUG5) {
			LOGGER(DEBUG5, "restrict_to_maybe_sat_prefix",
					"restricting to time domain to prefix satisfying constraint "+logger::formatted_to_string(con)+".");
		}

		std::vector<time_interval> itv_vec = get_maybe_satisfying_subdomains(
				con, err, false);

		time_interval it_pref = time_interval::empty();
		// accept the interval if it starts at the beginning of the time domain
		if (itv_vec.size() > 0
				&& is_MEQ(itv_vec.front().lower(), tdom.lower())) {
			// use the lower bound of the dom to avoid numerical difficulties
			it_pref = time_interval(tdom.lower(), itv_vec.front().upper());
		}
		sat_prefix = plif::intersect(sat_prefix, it_pref);
	}

	if (!plif::empty(sat_prefix)) {
		LOGGER(DEBUG5, "restrict_to_maybe_sat_prefix",
				"restricting to time interval ["+to_string(sat_prefix.lower())+","+to_string(sat_prefix.upper())+"]");
	} else {
		LOGGER(DEBUG5, "restrict_to_maybe_sat_prefix", "constraint unsat in prefix");
	}

	my_splif.restrict_to_subdomain(sat_prefix);
}

template<typename scalar_type>
void spacetime_flowpipe<scalar_type>::restrict_to_definitely_sat_prefix(
		const math::lin_constraint_system<scalar_type>& cons,
		const error_type& final_err) {
	using namespace math::numeric;

	LOGGERSW(DEBUG4, "restrict_to_definitely_sat_prefix",
			"restricting space-time flowpipe to definitely satisfying prefix");

	const time_interval& tdom(my_splif.get_time_domain());
	time_interval sat_prefix(my_splif.get_time_domain());

	// iterate with decreasing error

	scalar_type factor(16),scale(factor*factor*factor*factor);
	/** Track constraints to be checked. Those satisfied on the whole domain don't need to be checked. */
	math::lin_constraint_system<scalar_type> check_cons(cons);

	// get the prefix interval on which a constraint is satisfied
	typename math::lin_constraint_system<scalar_type>::iterator it =
			check_cons.begin();

	bool itv_empty = false;
	bool itv_singular = false;
	while (scale >= 1 && !itv_empty && !itv_singular
			&& check_cons.begin() != check_cons.end()) {
		error_type err = scale * final_err;
		// remap to domain if necessary
		math::lin_expression<scalar_type> expr = it->get_l();
		expr.reorder(domain());
		math::lin_constraint<scalar_type> con(expr, it->get_sign());
		IFLOGGER(DEBUG5) {
			LOGGER(DEBUG5, __FUNCTION__,
					"restricting to time domain to prefix satisfying constraint "+logger::formatted_to_string(con)+".");
		}

		std::vector<time_interval> itv_vec_maybe =
				my_splif.get_maybe_satisfying_subdomains(con, err, false);
		std::vector<time_interval> itv_vec_def =
				my_splif.get_definitely_satisfying_subdomains(con, err);

		time_interval it_pref_def = time_interval::empty();
		// accept the interval if it starts at the beginning of the time domain
		if (itv_vec_def.size() > 0
				&& is_MEQ(itv_vec_def.front().lower(), tdom.lower())) {
			// use the lower bound of the dom to avoid numerical difficulties
			it_pref_def = time_interval(tdom.lower(),
					itv_vec_def.front().upper());
		}
		time_interval it_pref_maybe = time_interval::empty();
		// accept the interval if it starts at the beginning of the time domain
		if (itv_vec_maybe.size() > 0
				&& is_MEQ(itv_vec_maybe.front().lower(), tdom.lower())) {
			// use the lower bound of the dom to avoid numerical difficulties
			it_pref_maybe = time_interval(tdom.lower(),
					itv_vec_maybe.front().upper());
		}

		// if definitely always satisfied, remove constraint from check list
		if (!plif::empty(sat_prefix)
				&& math::maybe(is_GE(it_pref_def.upper(), tdom.upper()))) {
			it = check_cons.erase(it);
			LOGGER(DEBUG7, __FUNCTION__, "reduced checking "+to_string(check_cons.size())+" constraints");
		} else {
			// restrict search interval to maybe satisfying prefix and check the next constraint
			sat_prefix = plif::intersect(sat_prefix, it_pref_maybe);
			// on the last turn, restrict to definite part
			if (scale <= scalar_type(1)) {
				sat_prefix = plif::intersect(sat_prefix, it_pref_def);
			}
			++it;
		}

		itv_empty = plif::empty(sat_prefix);
		itv_singular = math::numeric::is_MEQ(sat_prefix.lower(),
				sat_prefix.upper());

		// If we checked all constraints, restart with tighter error bound
		if (it == check_cons.end()) {
			it = check_cons.begin();
			scale = scale/factor;
			LOGGER(DEBUG7, __FUNCTION__, "reducing scale to "+to_string(scale));
		}

		my_splif.restrict_to_subdomain(sat_prefix);
	}

	if (!plif::empty(sat_prefix)) {
		LOGGER(DEBUG5, __FUNCTION__,
				"restricting to time interval ["+to_string(sat_prefix.lower())+","+to_string(sat_prefix.upper())+"]");
	} else {
		LOGGER(DEBUG5, __FUNCTION__,
				"constraint unsat in prefix");
	}
}

template<typename scalar_type>
std::vector<typename spacetime_flowpipe<scalar_type>::time_interval> spacetime_flowpipe<
		scalar_type>::get_maybe_satisfying_subdomains(
		const math::lin_constraint<scalar_type>& con, const error_type& err, bool widen_by_error) {
	return my_splif.get_maybe_satisfying_subdomains(con, err, widen_by_error);
}

template<typename scalar_type>
std::vector<typename spacetime_flowpipe<scalar_type>::time_interval> spacetime_flowpipe<
		scalar_type>::get_maybe_satisfying_subdomains(
		const math::lin_constraint_system<scalar_type>& cons,
		const error_type& err, bool widen_by_error) {
	LOGGER(DEBUG4, "get_maybe_satisfying_subdomains",
			"computing satisfying time intervals");
	IFLOGGER(DEBUG5) {
		LOGGER(DEBUG5, "get_maybe_satisfying_subdomains",
				"constraints: "+logger::formatted_to_string(cons)+".");
	}

	// start with a single universe interval
	std::vector<time_interval> itv_vec(1, get_time_domain());

	for (typename math::lin_constraint_system<scalar_type>::const_iterator it =
			cons.begin(); it != cons.end() && !itv_vec.empty(); ++it) {
		// remap to domain if necessary
		math::lin_expression<scalar_type> expr = it->get_l();
		expr.reorder(domain());
		math::lin_constraint<scalar_type> con(expr, it->get_sign());

		std::vector<time_interval> con_itv_vec =
				get_maybe_satisfying_subdomains(con, err, widen_by_error);


		// intersect the intvs for this constraint with the intervals gotten so far
		itv_vec = plif::intersect_intervals(itv_vec, con_itv_vec);

		using plif::operator<<;
		LOGGER_OS(DEBUG7,__FUNCTION__) << "sat intervals for " << con << ": " << con_itv_vec << " intersected: " << itv_vec;
//		std::cout << "sat intervals:" << con_itv_vec << " resulting in "<< itv_vec << std::endl;
	}

	return itv_vec;
}

template<typename scalar_type>
const typename spacetime_flowpipe<scalar_type>::annotated_plif* spacetime_flowpipe<
		scalar_type>::get_evolution(const direction& d) {
	return my_splif.get_evolution(d);
}

template<typename scalar_type>
math::tribool spacetime_flowpipe<scalar_type>::decide_outer_contains(
		spacetime_flowpipe<scalar_type>& f, bool refine) const {
	// only check containment if *this is locked
	if (is_refinement_locked()) {
		LOGGER(DEBUG4, "decide_outer_contains",
				"deciding containment");
		// only check if the maps are the same and the bloating is contained
		bool check = true;
		if (my_map && f.my_map)
			check = check && (*my_map == *f.my_map);
		else {
			check = check && (!my_map && !f.my_map);
		}
		if (check && my_W) {
			check = check && (f.my_W && math::definitely(my_W->contains(f.my_W)));
		} else {
			check = check && (!f.my_W);
		}
		if (!check) {
			return math::indeterminate();
		}

		return my_splif.decide_outer_contains(f.my_splif, refine);
	} else {
		return math::indeterminate();
	}
}

template<typename scalar_type>
math::tribool spacetime_flowpipe<scalar_type>::decide_outer_contains(
		const spacetime_flowpipe<scalar_type>& f) const {
	spacetime_flowpipe<scalar_type>& nonconst_f = const_cast<spacetime_flowpipe<scalar_type>&>(f);
	return decide_outer_contains(nonconst_f, false);
}

template<typename scalar_type>
void spacetime_flowpipe<scalar_type>::set_map(const affine_map& M) {
	my_map = affine_map_ptr(new affine_map(M));
}

template<typename scalar_type>
void spacetime_flowpipe<scalar_type>::set_bloating(
		const support_function_provider::const_ptr& W) {
	my_W = W;
}

template<typename scalar_type>
const typename spacetime_flowpipe<scalar_type>::affine_map_ptr& spacetime_flowpipe<
		scalar_type>::get_map() const {
	return my_map;
}

template<typename scalar_type>
support_function_provider::const_ptr spacetime_flowpipe<scalar_type>::get_initial_states() const {
	return my_splif.get_affine_support_flowpipe_description().X0;
}

template<typename scalar_type>
const support_function_provider::const_ptr& spacetime_flowpipe<scalar_type>::get_bloating() const {
	return my_W;
}

template<typename scalar_type>
bool spacetime_flowpipe<scalar_type>::is_refinement_locked() const {
	return my_refinement_locked;
}

template<typename scalar_type>
void spacetime_flowpipe<scalar_type>::set_refinement_locked() {
	my_refinement_locked = true;
}

template<typename scalar_type>
void spacetime_flowpipe<scalar_type>::restrict_to_subdomain(
		const time_interval& intv) {
	if (intv.is_empty()) {
		throw std::runtime_error("restricting to empty subdomain");
	}
	my_splif.restrict_to_subdomain(intv);
}

template<typename scalar_type>
const typename spacetime_flowpipe<scalar_type>::annotated_plif& spacetime_flowpipe<
		scalar_type>::get_or_compute_evolution(const direction& d,
		const error_type& err) {
	// if locked, redirect to get_evolution
	if (is_refinement_locked()) {
		const spacetime_plif::annotated_plif* p = get_evolution(d);
		if (p) {
			return *p;
		} else {
			throw basic_exception(
					"trying to compute evolution on locked spacetime_flowpipe");
		}
	} else {
		return my_splif.get_or_compute_evolution(d, err);
	}
}

template<typename scalar_type>
spacetime_flowpipe<scalar_type>::spacetime_flowpipe(
		const flowpipe_description& asfd, const time_interval& tdom) :
		index_to_variable_id_map_provider(asfd.dyn.domain()), my_splif(asfd,
				tdom), my_refinement_locked(false) {
	my_splif.set_simplify_concave(simplify_concave);
	my_splif.set_simplify_convex(simplify_convex);
}

template<typename scalar_type>
spacetime_flowpipe<scalar_type>::spacetime_flowpipe(spacetime_plif splif,
		affine_map_ptr map, support_function_provider::const_ptr W,
		bool refinement_locked, math::lin_constraint_system<scalar_type> cons) :
		index_to_variable_id_map_provider(
				splif.get_affine_support_flowpipe_description().dyn.domain()), my_splif(
				splif), my_map(map), my_W(W), my_refinement_locked(refinement_locked), my_outer_constraints(cons) {

}


template<typename scalar_type>
const variable& spacetime_flowpipe<scalar_type>::get_time_variable() const {
	return my_splif.get_time_variable();
}

template<typename scalar_type>
spacetime_flowpipe<scalar_type>::~spacetime_flowpipe<scalar_type>() {
}

template<typename scalar_type>
continuous_set_predicate::ptr spacetime_flowpipe<scalar_type>::get_predicate() const {
	throw std::runtime_error("missing implementation");
	return continuous_set_predicate::ptr();
}

template<typename scalar_type>
spacetime_flowpipe<scalar_type>* spacetime_flowpipe<scalar_type>::create_universe() const {
	return new spacetime_flowpipe<scalar_type>(
			my_splif.get_affine_support_flowpipe_description(),
			my_splif.get_time_domain());
}

template<typename scalar_type>
spacetime_flowpipe<scalar_type>* spacetime_flowpipe<scalar_type>::create_empty() const {
	return new spacetime_flowpipe<scalar_type>(
			my_splif.get_affine_support_flowpipe_description(),
			time_interval::empty());
}

template<typename scalar_type>
spacetime_flowpipe<scalar_type>* spacetime_flowpipe<scalar_type>::clone() const {
	return new spacetime_flowpipe<scalar_type>(*this);
}

template<typename scalar_type>
int spacetime_flowpipe<scalar_type>::get_memory() const {
	throw std::runtime_error("missing implementation");
	return 0;
}

template<typename scalar_type>
unsigned int spacetime_flowpipe<scalar_type>::get_dim() const {
	throw std::runtime_error("missing implementation");
	return 0;
}

template<typename scalar_type>
math::tribool spacetime_flowpipe<scalar_type>::is_empty() const {
	//throw std::runtime_error("missing implementation");
	if (my_splif.get_time_domain().is_empty()) {
		return true;
	} else {
		return math::indeterminate();
	}
}

template<typename scalar_type>
void spacetime_flowpipe<scalar_type>::embed_variables(
		const variable_id_set& id_set) {
	throw std::runtime_error("missing implementation");
}

template<typename scalar_type>
void spacetime_flowpipe<scalar_type>::existentially_quantify_variables(
		const variable_id_set& id_set) {
	throw std::runtime_error("missing implementation");
}

template<typename scalar_type>
void spacetime_flowpipe<scalar_type>::accept(
		dispatching::dispatcher<continuous_set_typelist>& d) const {
	d.dispatch(this);
//	throw std::runtime_error("missing implementation");
}

template<typename scalar_type>
void spacetime_flowpipe<scalar_type>::print(std::ostream& os) const {
	os << "spacetime_flowpipe(";
	os << "domain: " << domain();
 	os << ", locked: " << std::boolalpha << my_refinement_locked;
	os << ", samples: "; my_splif.print(os);
	if (my_map) {
		os << ", mapped by: " << *my_map;
	}
	if (my_W) {
		os << ", bloated by: " << *my_W;
	}
	os << ", outer constraints: " << my_outer_constraints;
	os << ")";
}

template<typename scalar_type>
spacetime_flowpipe<scalar_type>::spacetime_flowpipe() {

}

template<typename scalar_type>
support_function_provider::const_ptr spacetime_flowpipe<scalar_type>::get_last_states() const {
	return my_splif.get_last_states();
}

template<typename scalar_type>
const typename spacetime_flowpipe<scalar_type>::flowpipe_description& spacetime_flowpipe<scalar_type>::get_flowpipe_description() const {
	return my_splif.get_affine_support_flowpipe_description();
}



}

#endif /* SPACETIME_FLOWPIPE_HPP_ */
