/*
 * spacetime_scenario.h
 *
 *  Created on: Nov 3, 2012
 *      Author: notroot
 */

#ifndef SPACETIME_SCENARIO_H_
#define SPACETIME_SCENARIO_H_

#include <ostream>

#include "math/numeric/boost_interval_utility.h"
#include "core/scenarios/reachability_scenario.h"

namespace hybrid_automata {

reachability_scenario get_spacetime_scenario();

}


#endif /* SPACETIME_SCENARIO_H_ */
