/*
 * scenario_chooser.h
 *
 *  Created on: Jun 21, 2009
 *      Author: frehse
 */

#ifndef SCENARIO_CHOOSER_H_
#define SCENARIO_CHOOSER_H_

#include <string>

#include "core/scenarios/reachability_scenario.h"

namespace hybrid_automata {

class scenario_chooser {
public:
	static void set_scenario(const std::string& opt);
	static void set_scenario(reachability_scenario scen);
	static const reachability_scenario& get_scenario();
private:
	static reachability_scenario my_scenario;
	// forbid copying and assigning
	scenario_chooser(const scenario_chooser&);
	scenario_chooser& operator =(const scenario_chooser&);
};

}

#endif /* SCENARIO_CHOOSER_H_ */
