/*
 * hybrid_automaton_utility.cpp
 *
 *  Created on: Sep 21, 2009
 *      Author: frehse
 */

#include "core/hybrid_automata/hybrid_automaton_utility.h"

#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/location.h"
#include "core/hybrid_automata/automaton_cache.h"
#include "core/hybrid_automata/adapt_automaton_visitor.h"
#include "core/discrete/discrete_set_stl_set.h"
#include "core/discrete/singleton_set.h"
#include "core/symbolic_states/symbolic_state.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/continuous/continuous_set_operators.h"

namespace hybrid_automata {

/** Old version assuming that automata might not return a fully canonic version of their constraints */
location_constraint_set canonicalize_location_constraints_old(
		const location_constraint_set& lcons) {
	location_constraint_set newcons;
	for (location_constraint_set::const_iterator it = lcons.begin(); it
			!= lcons.end(); ++it) {
		location_constraint_set tmpcons;
		hybrid_automaton::ptr aut = hybrid_automaton_cache::get_automaton(
				it->first);
		if (!aut) {
			hybrid_automaton_cache::print(std::cerr);
			throw std::runtime_error("could not find automaton with id "
					+ int2string(it->first) + ".");
		}
		bool changed = aut->canonicalize_location_constraint(it->first,
				it->second, tmpcons);
		if (changed) {
			tmpcons = canonicalize_location_constraints_old(tmpcons);
		} else
			newcons.intersection_assign(tmpcons);
	}
	return newcons;
}
;

location_constraint_set canonicalize_location_constraints(
		const location_constraint_set& lcons) {
	location_constraint_set newcons;
	for (location_constraint_set::const_iterator it = lcons.begin(); it
			!= lcons.end(); ++it) {
		hybrid_automaton::ptr aut = hybrid_automaton_cache::get_automaton(
				it->first);
		if (!aut) {
			hybrid_automaton_cache::print(std::cerr);
			throw std::runtime_error("could not find automaton with id "
					+ int2string(it->first) + ".");
		}
		aut->canonicalize_location_constraint(it->first, it->second, newcons);
	}
	return newcons;
}

class canocalize_location_constraint_visitor: public adapt_discrete_set_visitor {
	virtual void dispatch(const discrete::discrete_set_stl_set* c) {
		discrete::discrete_set_stl_set::ptr d =
				discrete::discrete_set_stl_set::ptr(
						new discrete::discrete_set_stl_set());
		const discrete::discrete_set_stl_set::container_type& cont =
				c->get_stl_set();
		for (discrete::discrete_set_stl_set::container_type::const_iterator it =
				cont.begin(); it != cont.end(); ++it) {
			d->add(canonicalize_location_constraints(*it));
		}
		my_discrete_set = d;
		success = true;
	}
	virtual void dispatch(const discrete::singleton_set* c) {
		discrete::singleton_set::ptr d = discrete::singleton_set::ptr(
				new discrete::singleton_set(canonicalize_location_constraints(
						c->get_object())));
		my_discrete_set = d;
		success = true;
	}
};

bool canonicalize_location_constraints(symbolic_state_collection& s) {
	canocalize_location_constraint_visitor v;
	s.accept(v);
	return v.get_success();
}

bool adapt(adapt_automaton_visitor& v, symbolic_state_collection& sstates) {
	adapt_continuous_set_visitor_ptr cv = v.get_adapt_continuous_set_visitor();
	cv->reset();
	sstates.accept(*cv);
	return cv->get_success();
}

symbolic_state_collection_ptr intersect_with_invariants(
		const symbolic_state_collection& sstates,
		const hybrid_automaton& aut) {
	// for the result, create a collection of the same type as sstates
	symbolic_state_collection_ptr res_states = symbolic_state_collection::ptr(
			sstates.create());

	// iterate over symbolic states in collection
	for (symbolic_state_collection::const_iterator ssit = sstates.begin(); ssit
			!= sstates.end(); ++ssit) {
		const discrete::discrete_set::const_ptr& orig_dset = (*ssit)->get_discrete_set();
		const continuous::continuous_set::const_ptr& cset=(*ssit)->get_continuous_set();

//		std::cout << std::endl << "orig dset: " << orig_dset << std::endl;

		// for each symbolic state,
		// obtain discrete sets that are equivalent with respect to time elase
		typedef std::list<discrete::discrete_set_ptr> dset_list_type;
		dset_list_type dset_list = aut.get_time_equiv_locations(
				orig_dset, cset);
		for (dset_list_type::const_iterator dit = dset_list.begin(); dit
				!= dset_list.end(); ++dit) {
			// for each of the discrete sets,
			// construct a new symbolic state whose continuous set is intersected with the invariant
			// and add it to the result
			const discrete::discrete_set_ptr& dset = *dit;

			// Each element of dset corresponds to a set constraints,
			// each of which corresponds to a set of of locations.
			// Since they are all equivalent, we can just take the first location in the first set.
			location_id_set locs = aut.get_locations(*dset->begin());

// std::cout << std::endl << "dset: " << dset << ", locs: " << locs << std::endl;

			if (!locs.empty()) {
				// All locs in dset are supposed to have the same dynamics, so let's just choose the first one
				// and get its invariant
				location_id some_loc = *locs.begin();
				const continuous::continuous_set::const_ptr
						& inv =
								aut.get_location(some_loc)->get_time_constraints().get_invariant();

				// intersect the continuous part of the symbolic state
				// with the invariant
				continuous::continuous_set::ptr new_cset = continuous::compute_intersection(cset, inv);

				// add the corresponding symbolic state to the result set
				if (!math::definitely(new_cset->is_empty())) {
					res_states->add(symbolic_state::ptr(new symbolic_state(
							dset, new_cset)));
				}
			}
		}
	}

//std::cout << "before intersection:" << sstates << std::endl;
//std::cout << "after intersection:" << res_states << std::endl;
	return res_states;
}

location_constraint_set canonicalize_location_constraint(
		const hybrid_automaton::const_ptr& aut, const location_id& loc_id) {
	location_constraint_set lcons;

	location_constraint con(loc_id);
	aut->canonicalize_location_constraint(aut->get_id(), con, lcons);

	return lcons;
}

}
