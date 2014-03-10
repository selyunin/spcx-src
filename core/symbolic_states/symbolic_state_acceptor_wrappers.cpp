/*
 * symbolic_state_acceptor_wrappers.cpp
 *
 *  Created on: Sep 26, 2009
 *      Author: frehse
 */

#include "core/symbolic_states/symbolic_state_acceptor_wrappers.h"

//#include "../../abstract_framework/hybrid_automata/hybrid_automaton.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/pwl/passed_and_waiting_list.h"

namespace hybrid_automata {

plwl_sink::plwl_sink(passed_and_waiting_list& plwl, hybrid_automaton_ptr& aut) :
	my_plwl(plwl), my_aut(aut) {
}

void plwl_sink::accept(const symbolic_state_ptr& sstate) {
	my_plwl.add_passed(sstate, my_aut);
	my_plwl.add_waiting(sstate, my_aut);
}

symbolic_state_collection_sink::symbolic_state_collection_sink(symbolic_state_collection& scoll) :
	my_scoll(scoll) {
}

void symbolic_state_collection_sink::accept(const symbolic_state_ptr& sstate) {
	my_scoll.add(sstate);
}
}
