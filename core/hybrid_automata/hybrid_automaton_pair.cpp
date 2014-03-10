/*
 * hybrid_automaton_pair.cpp
 *
 *  Created on: Aug 26, 2009
 *      Author: frehse
 */

#include "core/hybrid_automata/hybrid_automaton_pair.h"
#include "core/hybrid_automata/composition_operators.h"
#include "core/discrete/discrete_set.h"
#include "core/symbolic_states/symbolic_state.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/discrete/singleton_set.h"

/** Forward declarations */
namespace continuous {
class continuous_set;
typedef boost::shared_ptr<const continuous_set> continuous_set_const_ptr;
}

namespace hybrid_automata {

using namespace discrete;
using namespace continuous;

bool hybrid_automaton_pair::find_location_id(location_id& loc, const location_id& id1,
		const location_id& id2) const {
	//std::cout << my_pair_to_loc_bimap.left << "!" << std::flush;
	pair_to_loc_bimap::left_map::const_iterator it = my_pair_to_loc_bimap.left.find(std::make_pair(
			id1, id2));
	if (it != my_pair_to_loc_bimap.left.end()) {
		loc = it->second;
		return true;
	} else
		return false;
}

location_id hybrid_automaton_pair::get_or_add_location_id(const location_id& id1,
		const location_id& id2) {
	assert(left_child);
	assert(right_child);

	location_id loc;
	if (!find_location_id(loc, id1, id2)) {
		// loc doesn't yet exist

		// add the appropriate location right now
		loc = add_location(compose_location(left_child->get_location(id1),
				right_child->get_location(id2)));

		// add the location to the bimap
		my_pair_to_loc_bimap.insert(pair_to_loc_bimap::value_type(std::make_pair(id1, id2), loc));
	}
	return loc;
}
// GF: This seems obsolete since location_constraint_sets handle locations in the
//     composition just fine
//
// /** Add the composition of the two symbolic states to res.
// *
// * @note res is passed as a parameter so that its implementation type is
// * determined by the caller.
// */
//void hybrid_automaton_pair::add_composed_states(symbolic_state_collection::ptr& res,
//		const symbolic_state::const_ptr& s1, const symbolic_state::const_ptr& s2) {
//
//	discrete_set::const_ptr d1 = s1->get_discrete_set();
//	discrete_set::const_ptr d2 = s2->get_discrete_set();
//
//	// The continuous set is the intersection of the two continuous_sets
//	continuous_set::ptr c =
//			compute_intersection(s1->get_continuous_set(), s2->get_continuous_set());
//
//	// for the locations build the cross product:
//	// enumerate all locations of s1 with all of s2
//	for (discrete_set::const_iterator it = d1->begin(); it != d1->end(); ++it) {
//		for (discrete_set::const_iterator jt = d2->begin(); jt != d2->end(); ++jt) {
//			// Intersect the location constraints
//			// Instantiate as a singleton discrete set
//			location_constraint_set lcs = *it;
//			lcs.intersection_assign(*jt);
//			if (!lcs.is_empty()) {
//				singleton_set::ptr d = singleton_set::ptr(new singleton_set(lcs));
//
//				symbolic_state::ptr new_s = symbolic_state::ptr(new symbolic_state(d, c));
//				res->add(new_s);
//
//				// Instantiating the locations in the set
//				// is not necessary
//				//location_id_set lis1=left_child->get_locations(*it);
//				//location_id_set lis2=right_child->get_locations(*it);
//			}
//		}
//	}
//}
//
//symbolic_state_collection::ptr hybrid_automaton_pair::compose(
//		const symbolic_state_collection::const_ptr& s1,
//		const symbolic_state_collection::const_ptr& s2) {
//	// create a collection of the same type as that of s1
//	symbolic_state_collection::ptr res = s1->create();
//	symbolic_state_collection::ptr new_states = s1->create();
//
//	// build the cross product: enumerate all symbolic states of s1 with all of s2
//	for (symbolic_state_collection::const_iterator it = s1->begin(); it != s1->end(); ++it) {
//		for (symbolic_state_collection::const_iterator jt = s2->begin(); jt != s2->end(); ++jt) {
//			new_states->clear();
//			add_composed_states(new_states, *it, *jt);
//			res->union_assign(new_states);
//		}
//	}
//	return res;
//}

void hybrid_automaton_pair::update_outgoing_transitions(const location_id& l) const {
	// check if location is already updates
	if (my_loc_outgoing_uptodate_map.find(l) == my_loc_outgoing_uptodate_map.end()) {
		// need to const_cast *this in order to update implementation
		hybrid_automaton_pair* nonconst_this = const_cast<hybrid_automaton_pair*> (this);

		//std::cout << "looking for " << l << std::endl << std::flush;
		loc_id_pair l1l2 = my_pair_to_loc_bimap.right.at(l);
		location_id l1 = l1l2.first;
		location_id l2 = l1l2.second;

		label_id_set labs = get_labels();
//std::cout << "labels:" << labs << " silent:" << named_label::silent_id() << std::endl;

		// Update for all labels
		for (label_id_set::const_iterator label_it = labs.begin(); label_it
				!= labs.end(); ++label_it) {
			const label_id& a = *label_it;

			// Check which automata are participating
			bool left_part = left_child->get_labels().find(a) != left_child->get_labels().end();
			bool right_part = right_child->get_labels().find(a) != right_child->get_labels().end();

			// note: the silent label is in both alphabets, so if it's the silent label,
			//       left_part and right_part are true.
			if (left_part || a == named_label::silent_id()) { // left_child is participating
				std::pair<transition_const_iterator, transition_const_iterator> out1 =
						left_child->get_outgoing_transitions(l1, a);

				// Cycle through outgoing transitions of left_child
				for (transition_const_iterator left_it = out1.first; left_it != out1.second; ++left_it) {
					transition::ptr t1 = left_child->get_transition(*left_it);
					location_id k1 = t1->get_target();

					if (right_part && a != named_label::silent_id()) { // right_child is participating
						std::pair<transition_const_iterator, transition_const_iterator> out2 =
								right_child->get_outgoing_transitions(l2, a);

						// Cycle through outgoing transitions of right_child
						for (transition_const_iterator right_it = out2.first; right_it
								!= out2.second; ++right_it) {
							transition::ptr t2 = right_child->get_transition(*right_it);
							location_id k2 = t2->get_target();

							jump_constraints tcon = compose_jump_constraints(
									t1->get_jump_constraints(), t2->get_jump_constraints());

							location_id target_loc = nonconst_this->get_or_add_location_id(k1, k2);
							transition::ptr new_trans(t1->create(l, a, tcon,
									target_loc));
							my_impl->add_transition(new_trans);
						}
					} else { // it's only the left child, right child remains in loc l2
						jump_constraints tcon = compose_autonomous_jump_constraints(
								t1->get_jump_constraints(), left_child, right_child);
						transition::ptr new_trans(t1->create(l, a, tcon,
								nonconst_this->get_or_add_location_id(k1, l2)));
						my_impl->add_transition(new_trans);
					}
				} // end for
			} // end if
			if ((!left_part && right_part) || a == named_label::silent_id()) {
				// it's only the right child, left child remains in loc l1
				std::pair<transition_const_iterator, transition_const_iterator> out2 =
						right_child->get_outgoing_transitions(l2, a);

				// Cycle through outgoing transitions of right_child
				for (transition_const_iterator right_it = out2.first; right_it != out2.second; ++right_it) {
					transition::ptr t2 = right_child->get_transition(*right_it);
					location_id k2 = t2->get_target();
					jump_constraints tcon = compose_autonomous_jump_constraints(
							t2->get_jump_constraints(), right_child, left_child);
					transition::ptr new_trans(t2->create(l, a, tcon,
							nonconst_this->get_or_add_location_id(l1, k2)));
					my_impl->add_transition(new_trans);
				} // end for
			} // if ... participating

		}

		// mark location as updated
		nonconst_this->my_loc_outgoing_uptodate_map[l] = true;
	} // if uptodate
	// otherwise don't do anything
}

hybrid_automaton_pair::hybrid_automaton_pair() :
	left_child(hybrid_automaton::ptr()), right_child(hybrid_automaton::ptr()) {
}

void hybrid_automaton_pair::init(hybrid_automaton::ptr hcomp, hybrid_automaton::ptr h1,
		hybrid_automaton::ptr h2) {
	assert(hcomp);
	assert(h1);
	assert(h2);

	left_child = h1;
	right_child = h2;

	hybrid_automaton_wrapper::set_impl(hcomp);

	// Note: The name of my_impl is irrelevant as it is never accessed.
	set_name(left_child->get_name() + "~" + right_child->get_name());
	// hcomp is the current (possibly incomplete) on-the-fly subset of the composition
	hcomp->set_name(left_child->get_name() + "~" + right_child->get_name() + "_OTF_SUBSET");

	/** Define variables. */
	variable_id_set vars = h1->get_variable_ids();
	variable_id_set var2 = h2->get_variable_ids();
	vars.insert(var2.begin(), var2.end());
	variable_id_set inp_vars = compute_input_variables(h1, h2);
	variable_id_set const_vars = compute_const_variables(h1, h2);
	add_variables(vars, inp_vars, const_vars);

	/** Define labels. */
	add_labels(h1->get_labels());
	add_labels(h2->get_labels());

	/** Define initial states as null pointer, will be instantiated when needed. */
	set_initial_states(symbolic_state_collection_ptr());
}

hybrid_automaton_pair::~hybrid_automaton_pair() {
}

const symbolic_state_collection_ptr& hybrid_automaton_pair::get_initial_states() const {
	if (!my_impl->get_initial_states()) {
		// not yet created
		symbolic_state_collection_ptr new_states = compose_symbolic_states(
				left_child->get_initial_states(), right_child->get_initial_states());
		// need to const_cast *this in order to update implementation
		//hybrid_automaton_pair* nonconst_this = const_cast<hybrid_automaton_pair*> (this);
		my_impl->set_initial_states(new_states);
	}
	return my_impl->get_initial_states();
}

void hybrid_automaton_pair::set_initial_states(const symbolic_state_collection::ptr& sstate_set) {
	//	// instantiate the states by obtaining the location_ids for it
	//	for (symbolic_state_collection::const_iterator it = sstate_set->begin(); it
	//			!= sstate_set->end(); ++it) {
	//		const discrete_set::const_ptr d = (*it)->get_discrete_set();
	//		for (discrete_set::const_iterator jt = d->begin(); jt != d->end(); ++jt) {
	//			get_locations(*jt); // this forces the instantiation for all locations in *jt
	//		}
	//	}

	my_impl->set_initial_states(sstate_set);
}

std::pair<hybrid_automaton_pair::transition_const_iterator,
		hybrid_automaton_pair::transition_const_iterator> hybrid_automaton_pair::get_outgoing_transitions(
		location_id l, label_id a) const {
	update_outgoing_transitions(l);
	return hybrid_automaton_wrapper::get_outgoing_transitions(l, a);
}

std::pair<hybrid_automaton_pair::location_const_iterator,
		hybrid_automaton_pair::location_const_iterator> hybrid_automaton_pair::get_locations() const {
	// Instantiate all locations of the product space
	std::pair<location_const_iterator, location_const_iterator> locs1 = left_child->get_locations();
	std::pair<location_const_iterator, location_const_iterator> locs2 =
			right_child->get_locations();

	// need to const_cast *this in order to update implementation
	hybrid_automaton_pair* nonconst_this = const_cast<hybrid_automaton_pair*> (this);

	for (location_const_iterator it = locs1.first; it != locs1.second; ++it) {
		for (location_const_iterator jt = locs2.first; jt != locs2.second; ++jt) {
			nonconst_this->get_or_add_location_id(*it, *jt);
		}
	}

	// Now we're sure every location exists; it suffices to forward
	// the call to the implementation
	return hybrid_automaton_wrapper::get_locations();
}

location_id_set hybrid_automaton_pair::get_locations(const location_constraint_set& lcons) const {
	location_id_set locset, set1, set2;
	location_id_set neglocs;

	// Add locations that might be defined over *this
	// shortcut if it's a single constraint that defines the top implementation
	location_constraint_set::const_iterator_pair itp = lcons.get_constraints(get_id());
	if (itp.first == itp.second) {
		itp = lcons.get_constraints(my_impl->get_id());
	}
	if (itp.first != itp.second) {
		const location_constraint& con = itp.first->second;
		if (con.get_sign()) {
			// exactly one location corresponds to the constraint.
			// since the id is already out there, it has to be known
			locset.insert(con.get_id());
			return locset;
		} else {
			// all locations except con->get_id() need to be added
			// we define therefore a negative set
			for (location_constraint_set::const_iterator it = itp.first; it != itp.second; ++it) {
				neglocs.insert(it->second.get_id());
			}

			// get all locations and subtract
			location_constraint_set emptycons;
			locset = get_locations(emptycons);
			if (!neglocs.empty()) {
				set_difference_assign(locset, neglocs);
			}
		}
	}

	// Add the product of locations of the children
	set1 = left_child->get_locations(lcons);
	set2 = right_child->get_locations(lcons);

	// need to const_cast *this in order to update implementation
	hybrid_automaton_pair* nonconst_this = const_cast<hybrid_automaton_pair*> (this);

	for (location_id_set::const_iterator it = set1.begin(); it != set1.end(); ++it) {
		for (location_id_set::const_iterator jt = set2.begin(); jt != set2.end(); ++jt) {
			locset.insert(nonconst_this->get_or_add_location_id(*it, *jt));
		}
	}

	return locset;
}

const hybrid_automaton_pair::loc_id_pair& hybrid_automaton_pair::get_child_locations(const location_id& l) const {
	pair_to_loc_bimap::right_const_iterator right_iter = my_pair_to_loc_bimap.right.find(l);
	if (right_iter == my_pair_to_loc_bimap.right.end()) {
		throw std::runtime_error("Could not find location_id " + to_string(l)
				+ " in get_child_locations(l)");
	}
	return right_iter->second;
}

bool hybrid_automaton_pair::canonicalize_location_constraint(const automaton_id& aut_id, const location_constraint& con,
		location_constraint_set& lcons) const {
	bool changed=false;

	// This the constraint is about *this, then look up child location pairs and propagate
	// Translate the location id to children's location ids
	if (aut_id == get_id() || aut_id == my_impl->get_id()) {
		loc_id_pair lip = get_child_locations(con.get_id());
		location_constraint_set childcons1;
		changed = changed
				|| left_child->canonicalize_location_constraint(
						left_child->get_id(),
						location_constraint(lip.first, con.get_sign()), lcons);
		changed = changed
				|| right_child->canonicalize_location_constraint(
						right_child->get_id(),
						location_constraint(lip.second, con.get_sign()), lcons);
	} else {
		// It could still be a child constraint. Try both.
		changed = changed
				|| left_child->canonicalize_location_constraint(
						aut_id, con, lcons);
		changed = changed
				|| right_child->canonicalize_location_constraint(
						aut_id, con, lcons);
	}
	return changed;
}

void hybrid_automaton_pair::accept(hybrid_automaton_visitor& v) {
	if (v.get_network_visits()
			== hybrid_automaton_visitor::all_components_and_subsets
			|| v.get_network_visits()
					== hybrid_automaton_visitor::all_base_components) {
		left_child->accept(v);
		right_child->accept(v);
	}
	if (v.get_network_visits()
			== hybrid_automaton_visitor::all_components_and_subsets
			|| v.get_network_visits()
					== hybrid_automaton_visitor::only_composition)
		hybrid_automaton_wrapper::accept(v);
}


}
