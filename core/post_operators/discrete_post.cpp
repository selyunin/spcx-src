/*
 * discrete_post.cpp
 *
 *  Created on: Aug 19, 2009
 *      Author: frehse
 */

#include "core/post_operators/discrete_post.h"

//#include "../continuous/continuous_set.h"
//#include "../discrete/discrete_set.h"
#include "core/discrete/singleton_set.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/hybrid_automaton_utility.h"
#include "core/hybrid_automata/location.h"
#include "core/hybrid_automata/location_constraint_set.h"
#include "core/hybrid_automata/location_constraint_set_operators.h"
#include "core/hybrid_automata/transition.h"
#include "core/symbolic_states/symbolic_state.h"
#include "core/symbolic_states/symbolic_state_collection.h"

namespace hybrid_automata {

void discrete_post::add_post_states_trans(
		const hybrid_automaton::const_ptr& aut,
		symbolic_state_collection_ptr& passed_result_set,
		symbolic_state_collection_ptr& waiting_result_set, const transition_id& trans,
		const symbolic_state::ptr& sstate) const {
	// Get the target location of the transition.
	if (sstate) {
		transition::ptr trans_ptr = aut->get_transition(trans);
		const location_id& s_id = trans_ptr->get_source();
		const location_id& t_id = trans_ptr->get_target();
		try {
			location_constraint_set source_cons = canonicalize_location_constraint(aut, s_id);
			location_constraint_set target_cons = canonicalize_location_constraint(aut, t_id);
			location_constraint_set source_cons_compact = compute_changes(target_cons,source_cons);
			location_constraint_set target_cons_compact = compute_changes(source_cons,target_cons);
			LOGGER_OS(HIGH,"discrete_post::add_post_states_trans") << "discrete post with label \"" <<  named_label::get_name(trans_ptr->get_label()) << "\" from " << source_cons_compact << " to " << target_cons_compact << std::endl << "source location: " << source_cons;

			location_constraint_set lcs = location_constraint_set(
					aut->get_id(), t_id);
			discrete::discrete_set::ptr dset = discrete::discrete_set::ptr(
					new discrete::singleton_set(lcs));
			// Get the invariant of the source location
			continuous::continuous_set::const_ptr s_inv = aut->get_location(
					s_id)->get_time_constraints().get_invariant();
			// Get the invariant of the target location
			continuous::continuous_set::const_ptr t_inv = aut->get_location(
					t_id)->get_time_constraints().get_invariant();
			const jump_constraints& jump = trans_ptr->get_jump_constraints();

			bool has_guard = (jump.get_guard());
			bool has_empty_guard = false;
			if (has_guard) {
				has_empty_guard = jump.get_guard()->is_empty();
			}
			if (!has_empty_guard) {
				const continuous::continuous_set::const_ptr& source_cset =
						sstate->get_continuous_set();
				continuous::continuous_set_collection csets = post(jump, s_inv,
						t_inv, source_cset);
				for (continuous::continuous_set_collection::iterator it =
						csets.begin(); it != csets.end(); ++it) {
					continuous::continuous_set::ptr cset = *it;
					if (cset && !math::definitely(cset->is_empty()))
						waiting_result_set->add(
								symbolic_state::ptr(
										new symbolic_state(dset, cset)));
				}
			}
		} catch (std::exception& e) {
			std::stringstream s;
			logger::copyfmt_to(s);
			s << canonicalize_location_constraint(aut, s_id);
			s << " to location ";
			s << canonicalize_location_constraint(aut, t_id);
			throw basic_exception(
					"could not apply discrete post for transition with id "
							+ to_string(trans) + ", label "
							+ named_label::get_name(trans_ptr->get_label())
							+ ", from location " + s.str() + ".", e);
		}
	}
}

void discrete_post::add_post_states(const hybrid_automaton::const_ptr& aut,
		symbolic_state_collection_ptr& passed_result_set,
		symbolic_state_collection_ptr& waiting_result_set, const label_id& lab,
		const symbolic_state::ptr& sstate) const {
	if (sstate) {
		std::list<transition_id> out_trans;
		out_trans = aut->get_outgoing_transitions(sstate->get_discrete_set(),
				lab);
		// iterate over transitions
		for (std::list<transition_id>::const_iterator trans_it =
			out_trans.begin(); trans_it != out_trans.end(); ++trans_it) {
			//			symbolic_state::ptr result_state = post(aut, *trans_it, sstate);
			//			if (result_state && !result_state->is_empty())
			//				result_set->add(result_state);
			add_post_states_trans(aut, passed_result_set, waiting_result_set, *trans_it, sstate);
		}
	}
}

void discrete_post::add_post_states(const hybrid_automaton::const_ptr& aut,
		symbolic_state_collection_ptr& passed_result_set,
		symbolic_state_collection_ptr& waiting_result_set,
		const label_id_set& lab_set, const symbolic_state::ptr& sstate) const {
	if (sstate) {
		// iterate over all labels
		for (label_id_set::const_iterator lab_it = lab_set.begin(); lab_it
				!= lab_set.end(); ++lab_it) {
			add_post_states(aut, passed_result_set, waiting_result_set, *lab_it, sstate);
		}
	}
}

void discrete_post::add_post_states(const hybrid_automaton::const_ptr& aut,
		symbolic_state_collection_ptr& passed_result_set,
		symbolic_state_collection_ptr& waiting_result_set,
		const symbolic_state::ptr& sstate) const {
	if (sstate) {
		label_id_set lab_set = aut->get_labels();
		//lab_set.insert(named_label::silent_id());
		add_post_states(aut, passed_result_set, waiting_result_set, lab_set, sstate);
	}
}

}
