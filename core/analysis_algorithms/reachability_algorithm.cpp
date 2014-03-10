#include "core/analysis_algorithms/reachability_algorithm.h"

#include "utility/basic_warning.h"
#include "utility/logger.h"
#include "utility/logger_stopwatch.h"

#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/scenarios/reachability_scenario.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/pwl/passed_and_waiting_list.h"
#include "core/post_operators/discrete_post.h"
#include "core/post_operators/continuous_post.h"
#include "core/symbolic_states/symbolic_state.h"

namespace hybrid_automata {

using namespace continuous;

reachability_algorithm::reachability_algorithm(
		const reachability_scenario& scen) :
	my_scen(scen) {
}

void reachability_algorithm::apply_post_ops(const symbolic_state_ptr & psstate) {
	// The waiting states from post_c are passed directly
	// to post_d, without passing by the waiting list.
	// That way, post_c can return its argument
	// for the waiting list without any overhead (identity operator).

	//	std::cout << "s: " << psstate << std::endl;
	symbolic_state_collection::ptr passed_states = my_scen.create_symbolic_state_collection();
	symbolic_state_collection::ptr cwaiting_states =my_scen.create_symbolic_state_collection();

	my_post_c->add_post_states(my_aut, passed_states, cwaiting_states, psstate);
	my_plwl->add_passed(passed_states, my_aut);

	passed_states = my_scen.create_symbolic_state_collection();
	symbolic_state_collection::ptr dwaiting_states =
			my_scen.create_symbolic_state_collection();
	my_post_d->add_post_states(my_aut, passed_states, dwaiting_states, cwaiting_states);
	my_plwl->add_waiting(dwaiting_states, my_aut);
	my_plwl->add_passed(passed_states, my_aut);

}

bool reachability_algorithm::iter_check(int iter_count) {
	if (my_scen.get_iter_max() >= 0) {
		return iter_count < my_scen.get_iter_max();
	} else
		return true;
}

symbolic_state_collection::ptr reachability_algorithm::reach(
		hybrid_automaton::ptr H,
		const symbolic_state_collection::const_ptr& start_states) {
	LOGGERSW(LOW,"reachability_algorithm::reach", "Computing reachable states");

	my_plwl = my_scen.create_passed_and_waiting_list();
	my_aut = H;
	symbolic_state_ptr psstate;

	if (start_states->is_empty()) {
		basic_warning("reachability algorithm", "initial set empty",
				basic_warning::UNUSUAL_INPUT);
	}


	// waiting list is initialized with the initial states of the hybrid
	// automaton H.
	symbolic_state_collection::ptr ini_waiting_list = my_scen.create_symbolic_state_collection();

	ini_waiting_list->copy(start_states);

	my_plwl->set_waiting_list(ini_waiting_list);

	// get the post operators from the scenario
	my_post_c = my_scen.create_continuous_post_operator();
	my_post_c->set_time_horizon(my_scen.get_time_horizon());
	my_post_c->set_sampling_time(my_scen.get_sampling_time());
	my_post_d = my_scen.create_discrete_post_operator();

	// @todo: intersect the initial states with the invariants?

	passed_and_waiting_list::info pinf;

	int iter_count = 0;
	bool wl_empty = my_plwl->is_waiting_list_empty();

	while (iter_check(iter_count) && !wl_empty) {
		LOGGERSWNC(MEDIUM, "reachability_algorithm::reach", "Iteration " + to_string(iter_count));
		logger::logger_id sw_id = logger::get_last_id();

		// Pick a state from the waiting list
		psstate = my_plwl->pick();
		assert(psstate); // there should be an actual state

		// Apply the post operators to the state
		apply_post_ops(psstate);
		++iter_count;

		// Check if waiting list is empty
		wl_empty = my_plwl->is_waiting_list_empty();

		// get and display information about the progress
		pinf = my_plwl->get_info();
		LOGGER_ATTACH(MEDIUM, "reachability_algorithm::reach", to_string(
				pinf.pl_size) + " sym states passed, "
				+ to_string(pinf.wl_size) + " waiting ", sw_id);
	}
	if (!wl_empty) {
		LOGGER(LOW, "reachability_algorithm::reach","Performed max. number of iterations (" + to_string(iter_count) + ") without finding fixpoint.");
	} else {
		LOGGER(LOW, "reachability_algorithm::reach", "Found fixpoint after " + to_string(iter_count) + " iterations.");
	}
	return my_plwl->get_passed_list();
}

symbolic_state_collection::ptr reachability_algorithm::reach(
		hybrid_automaton::ptr H) {
	if(!H)
		std::runtime_error("Reach called with null automaton.\n");
	if(!(H->get_initial_states()))
		std::runtime_error("Reach called with null initial states. \n");
	return reach(H, H->get_initial_states());
}

reachability_algorithm::~reachability_algorithm() {
}

}
