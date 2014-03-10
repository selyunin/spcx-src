/*
 * continuous_post.cpp
 *
 *  Created on: Sep 1, 2009
 *      Author: frehse
 */

#include "core/post_operators/continuous_post.h"

#include "core/continuous/continuous_set.h"
#include "core/discrete/discrete_set.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/hybrid_automaton_utility.h"
#include "core/hybrid_automata/location.h"
#include "core/symbolic_states/symbolic_state.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/continuous/continuous_set_collection.h"

namespace hybrid_automata {

continuous::continuous_set::ptr continuous_post::post(
		const hybrid_automaton::const_ptr& aut,
		const discrete::discrete_set::const_ptr& dset,
		const continuous::continuous_set::const_ptr& cset) const {
	if (cset) {
		location_id_set locs = aut->get_locations(*dset->begin());
		if (locs.empty()) {
			return continuous::continuous_set::ptr();
		}
		// All locs in dset are supposed to have the same dynamics,
		// so let's just choose the dynamics of the first location
		// in the set
		const location_id& first_loc_id = *locs.begin();
		const location_ptr& first_loc = aut->get_location(first_loc_id);
		const time_constraints& tcons = first_loc->get_time_constraints();

		continuous::continuous_set::ptr res;

		try {
			LOGGER_OS(HIGH,"continuous_post::post") << "time elapse in " << canonicalize_location_constraint(aut,first_loc_id);

			res = post(tcons, cset);
		} catch (std::exception& e) {
			std::stringstream s;
			logger::copyfmt_to(s);
			s << canonicalize_location_constraint(aut, first_loc_id);
			throw basic_exception("Error computing time elapse in component "
					+ aut->get_name() + ", location " + s.str(), e);
		}
		return res;
	} else {
		return continuous::continuous_set::ptr();
	}
}

void continuous_post::add_post_states(const hybrid_automaton::const_ptr& aut,
		symbolic_state_collection_ptr& passed_result_set,
		symbolic_state_collection_ptr& waiting_result_set,
		const symbolic_state::ptr& sstate) const {
	typedef std::list<discrete::discrete_set::ptr> dlist_type;
	if (sstate) {
		dlist_type d_list = aut->get_time_equiv_locations(
				sstate->get_discrete_set(), sstate->get_continuous_set());
		for (dlist_type::const_iterator d_it = d_list.begin(); d_it
				!= d_list.end(); ++d_it) {
			continuous::continuous_set::ptr cset = post(aut, *d_it,
					sstate->get_continuous_set());
			// if collection, add individual elements
			if (continuous::continuous_set_collection::ptr csets =
					boost::dynamic_pointer_cast<
							continuous::continuous_set_collection>(cset)) {
				for (continuous::continuous_set_collection::iterator it =
						csets->begin(); it != csets->end(); ++it) {
					continuous::continuous_set::ptr set = *it;
					if (set && !math::definitely(set->is_empty())) {
						passed_result_set->add(
								symbolic_state::ptr(
										new symbolic_state(*d_it, set)));
						waiting_result_set->add(
								symbolic_state::ptr(
										new symbolic_state(*d_it, set)));
					}
				}
			} else if (cset) {
				passed_result_set->add(
						symbolic_state::ptr(new symbolic_state(*d_it, cset)));
				waiting_result_set->add(
						symbolic_state::ptr(new symbolic_state(*d_it, cset)));
			}
		}
	}
}

void continuous_post::set_time_horizon(double th) {
	my_time_horizon = th;
}

double continuous_post::get_time_horizon() const {
	return my_time_horizon;
}

void continuous_post::set_sampling_time(double ts) {
	my_sampling_time = ts;
}

double continuous_post::get_sampling_time() const {
	return my_sampling_time;
}
}
