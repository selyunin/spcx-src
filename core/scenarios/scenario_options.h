/*
 * scenario_options.h
 *
 *  Created on: Nov 12, 2009
 *      Author: frehse
 */

#ifndef SCENARIO_OPTIONS_H_
#define SCENARIO_OPTIONS_H_

#include "application/options.h"

namespace options {

void add_scenario_options();
bool check_scenario_options(options::options_processor::variables_map& vmap);
bool apply_scenario_options_wo_system(options::options_processor::variables_map& vmap);
bool apply_scenario_options_with_system(options::options_processor::variables_map& vmap);

}

#endif /* SCENARIO_OPTIONS_H_ */
