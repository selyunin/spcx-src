/*
 * reachability_options.h
 *
 *  Created on: Sep 28, 2009
 *      Author: frehse
 */

#ifndef REACHABILITY_OPTIONS_H_
#define REACHABILITY_OPTIONS_H_

#include "application/options.h"
#include "io/common_output/output_formatter.h"

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
class symbolic_state_collection;
typedef boost::shared_ptr<symbolic_state_collection> symbolic_state_collection_ptr;
typedef boost::shared_ptr<const symbolic_state_collection> symbolic_state_collection_const_ptr;
}

namespace options {

void execute_reachability(hybrid_automata::hybrid_automaton_ptr aut_net,
		hybrid_automata::symbolic_state_collection_ptr ini_states,
		options::options_processor::variables_map& vmap, std::vector<io::output_formatter*> of);

void add_reachability_options();
bool check_reachability_options(options::options_processor::variables_map& vmap);
bool apply_reachability_options(options::options_processor::variables_map& vmap);

}
#endif /* REACHABILITY_OPTIONS_H_ */
