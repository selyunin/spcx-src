#include "core/hybrid_automata/explicit_automaton.h"

//#include "hybrid_automaton_visitor.h" // to avoid circularities, it is only included here
#include "core/hybrid_automata/location.h"
#include "core/hybrid_automata/location_constraint_set.h"
#include "core/hybrid_automata/transition.h"

#include "utility/basic_exception.h"
#include "core/predicates/dot_context_lookup.h"

namespace hybrid_automata {

using namespace continuous;

explicit_automaton::explicit_automaton() {
}

explicit_automaton::explicit_automaton(std::string new_name) :
	hybrid_automaton(new_name) {
}

explicit_automaton::~explicit_automaton() {
}

explicit_automaton* explicit_automaton::create() const {
	return new explicit_automaton();
}

const symbolic_state_collection_ptr& explicit_automaton::get_initial_states() const {
	return my_initial_states;
}

/** The transition_const_iterator returned by get_outgoing_transitions() is implemented by an iterator
 * to loc_and_label_to_transition_id_map, but must dereference to a transition::ptr.
 * We overload the normal forward_iterator_wrapper in order to dereference as intended.
 */
std::pair<explicit_automaton::transition_const_iterator,
		explicit_automaton::transition_const_iterator> explicit_automaton::get_outgoing_transitions(
		location_id l, label_id a) const {
	typedef loc_and_label_to_transition_id_map::const_iterator MI;
	std::pair<MI, MI> r = my_outgoing_transitions.equal_range(std::make_pair(l, a));

	transition_const_iterator tb = transition_const_iterator::create<
			simple_iterators::second_in_pair_dereferencing<MI> >(r.first);
	transition_const_iterator te = transition_const_iterator::create<
			simple_iterators::second_in_pair_dereferencing<MI> >(r.second);
	return std::make_pair(tb, te);
}

void explicit_automaton::add_variable(const variable_id& vid, bool is_input, bool is_const) {
	my_variables.insert(vid);
	if (is_input) my_input_variables.insert(vid);
	if (is_const) my_const_variables.insert(vid);
}
const variable_id_set& explicit_automaton::get_variable_ids() const {
	return my_variables;
}
const variable_id_set& explicit_automaton::get_input_variables() const {
	return my_input_variables;
}
const variable_id_set& explicit_automaton::get_const_variables() const {
	return my_const_variables;
}

void explicit_automaton::add_label(const label_id& lab) {
//	if (lab == named_label::silent_id()) {
//		std::stringstream ss;
//		named_label::print_named_label_cache(ss);
//		throw basic_exception("silent label with id "+to_string(lab)+" should not be added to alphabet. Current global labels:\n"+ss.str());
//	}
	my_labels.insert(lab);
}

location_id explicit_automaton::add_location(const location::ptr& loc) {
	location_id lid = my_locations.insert(loc);
	my_location_name_to_id_map[loc->get_name()] = lid;
	return lid;
}

transition_id explicit_automaton::add_transition(const transition::ptr& trans, bool check_emptiness) {
	// check that the target location's invariant is not empty
	// and that the guard set is not empty

	// check that label is in alphabet
	label_id t_lab = trans->get_label();
	if (get_labels().find(t_lab)==get_labels().end())
		throw basic_exception("Internal error, adding transition with unkown label \""+named_label::get_name(t_lab)+"\"");

	transition_id tid = 0;
	bool
			inv_empty =
					check_emptiness
							&& get_location(trans->get_target())->get_time_constraints().get_invariant()
							&& math::definitely(get_location(trans->get_target())->get_time_constraints().get_invariant()->is_empty());
	bool guard_empty = check_emptiness
			&& trans->get_jump_constraints().get_guard()
			&& math::definitely(trans->get_jump_constraints().get_guard()->is_empty());
	if (!inv_empty) {
		if (!guard_empty) {
			tid = my_transitions.insert(trans);
			my_outgoing_transitions.insert(
					std::make_pair(
							std::make_pair(trans->get_source(),
									t_lab), tid));
		} else
			LOGGER(DEBUG7,__FUNCTION__,"skipping transition with label "+named_label::get_name(trans->get_label())+", guard empty");
	} else
		LOGGER(DEBUG7,__FUNCTION__,"skipping transition with label "+named_label::get_name(trans->get_label())+", target inv empty");

	return tid;
}

void explicit_automaton::set_initial_states(const symbolic_state_collection_ptr& sstate_set) {
	my_initial_states = sstate_set;
}

const label_id_set& explicit_automaton::get_labels() const {
	return my_labels;
}

transition::ptr explicit_automaton::get_transition(const transition_id& id) const {
	return my_transitions.get_ptr(id);
}
;

location::ptr explicit_automaton::get_location(const location_id& id) const {
	return my_locations.get_ptr(id);
}

location_id explicit_automaton::get_location_id(std::string loc_name) const {
	location_name_to_id_map::const_iterator it = my_location_name_to_id_map.find(loc_name);
	if (it == my_location_name_to_id_map.end())
		throw basic_exception("Location '" + loc_name + "' does not exist in automaton "+get_name()+".");
	return it->second;
}

std::pair<explicit_automaton::location_const_iterator, explicit_automaton::location_const_iterator> explicit_automaton::get_locations() const {
	//	forward_iterator_wrapper<const location::ptr,location_collection::const_iterator,second_in_pair_dereferencing> te(r.second);
	location_const_iterator tb = location_const_iterator::create<
			simple_iterators::first_in_pair_dereferencing<location_collection::const_iterator,
					location_const_iterator::value_type> >(my_locations.begin());
	location_const_iterator te = location_const_iterator::create<
			simple_iterators::first_in_pair_dereferencing<location_collection::const_iterator,
					location_const_iterator::value_type> >(my_locations.end());
	return std::make_pair(tb, te);
}

location_id_set explicit_automaton::get_locations(const location_constraint_set& lcons) const {
	location_id_set locs;
	// get the constraint
	location_constraint_set::const_iterator_pair itp = lcons.get_constraints(get_id());
	// if there is no constraint on this aut, add all locations
	if (itp.first == itp.second) {
		for (location_collection::const_iterator it = my_locations.begin(); it
				!= my_locations.end(); ++it) {
			locs.insert(it->first);
		}
	} else {
		for (location_constraint_set::const_iterator con_it = itp.first; con_it != itp.second; ++con_it) {
			const location_constraint& con = con_it->second;
			if (con.get_sign() == false) {
				for (location_collection::const_iterator it = my_locations.begin(); it
						!= my_locations.end(); ++it) {
					if (con.get_id() != it->first) {
						locs.insert(it->first);
					}
				}
			} else {
				locs.insert(con.get_id());
			}
		}
	}
	return locs;
}

bool explicit_automaton::canonicalize_location_constraint(const automaton_id& aut_id, const location_constraint& con,
		location_constraint_set& lcons) const {
	// No need to modify the constraint.
//std::cout << "aut "+get_name()+", id "+to_string(get_id())+", adding atomic constraint "+to_string(aut_id).

	// don't add useless constraints on automata with just one location
	if (aut_id == get_id() && my_locations.size() > 1 ) {
		lcons.add_constraint(aut_id,con);
	}
	return false;
}

//discrete::discrete_set::ptr explicit_automaton::post(const transition_id& trans,
//		const discrete::discrete_set::const_ptr& dset) {
//	// the transition is an explicit_transition
//	explicit_transition::ptr tr=get_explicit_transition(trans);
//	// create a discrete set of the same type as dset
//	discrete::discrete_set::ptr res_set=dset->create();
//
//	for (discrete::discrete_set::const_iterator it=dset->begin();it != dset->end();++it) {
//		if (tr->get_source()==*it) {
//			res_set->union_assign(tr->get_target());
//		}
//	}
//	return res_set;
//}

void explicit_automaton::print(std::ostream& os) const {

	print_visitor pv(os);

	explicit_automaton& const_this = const_cast<explicit_automaton&> (*this);
	const_this.accept(pv);

	/* Print initial states */
}

void explicit_automaton::accept(hybrid_automaton_visitor& v) {
	try {

		v.prologue(*this);

		//"outgoing_transitions" corresponds to the CIF format
		if(v.get_visiting_order() == hybrid_automaton_visitor::outgoing_transitions) {

			/* Visit locations */
			for (location_collection::const_iterator loc_it = my_locations.begin(); loc_it
					!= my_locations.end(); ++loc_it) {

				v.visit(*this, *loc_it->second, loc_it->first);
				/* Visit transitions */
				for (transition_collection::const_iterator t_it = my_transitions.begin(); t_it
					!= my_transitions.end(); ++t_it) {
					hybrid_automata::transition& t = *t_it->second;
					if(t.get_source() == loc_it->first)
							v.visit(*this, *t_it->second, t_it->first);
					}
				}
			}
		//else is for SX format
		else {
		/* Visit locations */
		for (location_collection::const_iterator it = my_locations.begin(); it
				!= my_locations.end(); ++it) {
			v.visit(*this, *it->second, it->first);
		}

		/* Visit transitions */
		for (transition_collection::const_iterator it = my_transitions.begin(); it
				!= my_transitions.end(); ++it) {
			v.visit(*this, *it->second, it->first);
			}
		}

		v.epilogue(*this);
	} catch (std::exception& e) {
		throw basic_exception("Failure in automaton " + dot_context::context_free_name(get_name()) + ".", e);
	}

}


}
