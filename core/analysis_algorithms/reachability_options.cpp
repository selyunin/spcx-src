/*
 * reachability_options.cpp
 *
 *  Created on: Sep 28, 2009
 *      Author: frehse
 */

#include "core/analysis_algorithms/reachability_options.h"

#include <iostream>
#include <boost/algorithm/string/trim.hpp>

#include "utility/stl_helper_functions.h"
#include "utility/logger_stopwatch.h"

//#include "../abstract_framework/hybrid_automata/automaton_cache.h"
#include "core/hybrid_automata/automaton_cache.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/hybrid_automaton_utility.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/scenarios/reachability_scenario.h"
#include "core/scenarios/scenario_chooser.h"
#include "core/scenarios/parse_with_scenario.h"
#include "core/analysis_algorithms/reachability_algorithm.h"
#include "core/pwl/passed_and_waiting_list.h"

#include "application/options.h"
#include "io/common_input/input_options.h" // get_system_name
#include "io/common_output/output_options.h"

namespace options {

void execute_reachability(hybrid_automata::hybrid_automaton::ptr aut_net,
		hybrid_automata::symbolic_state_collection::ptr ini_states,
		options::options_processor::variables_map& vmap, std::vector<io::output_formatter*> of) {
	const hybrid_automata::reachability_scenario& scen =
			hybrid_automata::scenario_chooser::get_scenario();

	// Define forbidden states if provided
	hybrid_automata::symbolic_state_collection::ptr bad_states =
			hybrid_automata::symbolic_state_collection::ptr();
	std::string bad_statestring;

	//std::cout << "forbidden states are:" << std::endl;

	if (options::options_processor::get_string_option(vmap,"forbidden",bad_statestring)) {
		if (bad_statestring != "") {
			std::string sys_name = get_system_name(vmap);

			parser::parse_policy ppol = parser::parse_policy::SX_policy();
			// accept unknown vars so that subcomponent variables can be specified too
			ppol.add_unknown_vars = true;

			try {
				bad_states
						= hybrid_automata::parse_symbolic_state_collection_with_scenario(
								scen, bad_statestring, sys_name,
								ppol);
			} catch (std::exception& e) {
				throw basic_exception("Failed to define forbidden states.", e);
			}

			hybrid_automata::canonicalize_location_constraints(*bad_states);
			/* DEBUG*/
			//std::cout << "bad_states:" << *bad_states << std::endl;
		}
	}

	//std::cout << "Checking " << sys_name << " for states " << bad_statestring << "." << std::flush;

	hybrid_automata::reachability_algorithm alg(scen);
	hybrid_automata::symbolic_state_collection::ptr states;
	states = alg.reach(aut_net, ini_states);

	//std::cout << "reachable:" << states << std::endl;
	if (states) {
		hybrid_automata::canonicalize_location_constraints(*states);
	}

	if (bad_states) {
//std::cout << "bad:" << bad_states << std::endl;
		bad_states->intersection_assign(states);
//std::cout << "bad&reachable:" << bad_states << std::endl;
		if (bad_states->is_empty())
			LOGGER(LOW, "execute_reachability", "Forbidden states are not reachable.");
		else
			LOGGER(LOW, "execute_reachability", "Forbidden states are reachable.");

		LOGGERSW(LOW, "execute_reachability", "Output of reachable bad states");
		for (unsigned int i = 0; i < of.size(); ++i) {
			of[i]->output(*bad_states);
		}
	} else {
		LOGGERSW(LOW, "execute_reachability", "Output of reachable states");
		for (unsigned int i = 0; i < of.size(); ++i) {
			of[i]->output(*states);
		}
	}
}

void add_reachability_options() {
	options::options_processor::config.add_options()(
			"forbidden",
			boost::program_options::value<std::string>(),
			"Defines a set of forbidden states. If provided, the states that intersect with the given forbidden states are output. Takes as argument the forbidden states.");
	options::options_processor::config.add_options()("iter-max",
			boost::program_options::value<std::string>(),
			"Set the maximum number of iterations for the reachability algorithm.");
	options::options_processor::config.add_options()(
			"PLWL",
			boost::program_options::value<std::string>(),
			"Options for the PLWL:\n- MERGEP_ON (OFF) : remove states on the passed list if they are contained in new states (on by default).");
}

bool check_reachability_options(options::options_processor::variables_map& vmap) {
	return true;
}
bool apply_reachability_options(options::options_processor::variables_map& vmap) {
	hybrid_automata::reachability_scenario scen =
			hybrid_automata::scenario_chooser::get_scenario();

	std::string s;
	if (options::options_processor::get_string_option(vmap,"PLWL",s)) {
		scen.get_passed_and_waiting_list()->set_option(s);
	}
	if (options::options_processor::get_string_option(vmap,"iter-max",s)) {
		int iter;
		try {
			iter = from_string<int> (s);
		} catch ( std::exception& e ) {
			throw basic_exception("Failed to parse option iter-max = \""+s+"\" as an integer.");
		}
		scen.set_iter_max(iter);
	}

	hybrid_automata::scenario_chooser::set_scenario(scen);
	return true;
}

}
