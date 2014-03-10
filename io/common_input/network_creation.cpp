/*
 * network_creation.cpp
 *
 *  Created on: Sep 7, 2009
 *      Author: frehse
 */

#include "io/common_input/network_creation.h"

#include "core/scenarios/scenario_chooser.h"

#include "core/hybrid_automata/automaton_cache.h"
#include "core/hybrid_automata/hybrid_automaton_network.h"

namespace parser {
namespace automaton_parser {

using namespace hybrid_automata;

hybrid_automata::hybrid_automaton_network_ptr create_automaton_network(const std::string& name) {
	const reachability_scenario& scen = scenario_chooser::get_scenario();
	hybrid_automaton_network_ptr aut_net = scen.create_hybrid_automaton_network();
	aut_net->set_name(name);
	hybrid_automata::hybrid_automaton_cache::add_automaton(aut_net);
	return aut_net;
}

hybrid_automata::hybrid_automaton::ptr create_automaton(const std::string& name) {
	const hybrid_automata::reachability_scenario& scen =
			hybrid_automata::scenario_chooser::get_scenario();
	hybrid_automata::hybrid_automaton::ptr aut = scen.create_hybrid_automaton();
	aut->set_name(name);
	hybrid_automata::hybrid_automaton_cache::add_automaton(aut);
	return aut;
}

}
}
