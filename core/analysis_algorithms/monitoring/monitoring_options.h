/*
 * monitoring_options.h
 *
 *  Created on: Sep 28, 2009
 *      Author: frehse
 */

#ifndef MONITORING_OPTIONS_H_
#define MONITORING_OPTIONS_H_

#include "application/options.h"
#include "io/common_output/output_options.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/symbolic_states/symbolic_state_collection.h"


namespace options {
void add_monitoring_options();
bool check_monitoring_options(options::options_processor::variables_map& vmap);
bool apply_monitoring_options(options::options_processor::variables_map& vmap);

/** Returns true if monitoring is to be carried out instead of normal reachability.
 */
bool use_monitoring(options::options_processor::variables_map& vmap);

/** Start monitoring process, parsing commands from std::in. */
void monitor(hybrid_automata::hybrid_automaton::ptr aut_net,
		hybrid_automata::symbolic_state_collection::ptr ini_states,
		options::options_processor::variables_map& vmap, std::vector<io::output_formatter*> of);
}
#endif /* MONITORING_OPTIONS_H_ */
