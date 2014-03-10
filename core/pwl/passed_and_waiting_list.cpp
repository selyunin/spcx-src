/*
 * passed_and_waiting_list.cpp
 *
 *  Created on: Sep 2, 2009
 *      Author: frehse
 */

#include "core/pwl/passed_and_waiting_list.h"

#include "core/symbolic_states/symbolic_state_collection.h"

namespace hybrid_automata {


void passed_and_waiting_list::add_passed(symbolic_state_collection_ptr sstate_set,
		const hybrid_automaton_ptr H) {
	for (symbolic_state_collection::const_iterator it = sstate_set->begin(); it
			!= sstate_set->end(); ++it) {
		add_passed(*it, H);
	}
}

void passed_and_waiting_list::add_waiting(symbolic_state_collection_ptr sstate_set,
		const hybrid_automaton_ptr H) {
	for (symbolic_state_collection::const_iterator it = sstate_set->begin(); it
			!= sstate_set->end(); ++it) {
		add_waiting(*it, H);
	}
}

void passed_and_waiting_list::print(std::ostream& os) const {
	os << "Passed: " << get_passed_list() << std::endl;
	os << "Waiting: " << get_waiting_list() << std::endl << std::flush;
}

void passed_and_waiting_list::set_option(std::string opt) {
	throw std::runtime_error(
			"passed_and_waiting_list: no default implementation of option "
					+ opt + ".");
}
}
