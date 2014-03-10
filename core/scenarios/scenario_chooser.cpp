/*
 * scenario_chooser.cpp
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#include "core/scenarios/scenario_chooser.h"

#include <stdexcept>
#include "core/scenarios/phaver_scenario.h"
#include "core/scenarios/spacetime_scenario.h"
#include "core/scenarios/support_fun_scenario.h"
#include "core/scenarios/simulation_scenario.h"


namespace hybrid_automata {

void scenario_chooser::set_scenario(const std::string& opt) {
	if (opt == "phaver") {
		my_scenario = get_phaver_scenario();
	} else 	if (opt == "stc") {
		my_scenario = get_spacetime_scenario();
	} else 	if (opt == "supp") {
		my_scenario = get_support_fun_scenario();
	} else if (opt =="simu"){
		my_scenario = get_simulation_scenario();
	}else{
		throw std::runtime_error("unknown scenario "+opt);
	}
}

void scenario_chooser::set_scenario(reachability_scenario scen) {
	my_scenario=scen;
}

const reachability_scenario& scenario_chooser::get_scenario() {
	return my_scenario;
}

reachability_scenario scenario_chooser::my_scenario = get_phaver_scenario();

}
