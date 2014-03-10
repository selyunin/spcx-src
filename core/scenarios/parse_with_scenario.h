/*
 * parse_with_scenario.h
 *
 *  Created on: Sep 23, 2009
 *      Author: frehse
 */

#ifndef PARSE_WITH_SCENARIO_H_
#define PARSE_WITH_SCENARIO_H_

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "io/common_input/parse_policy.h"

/** Forward declarations */
namespace hybrid_automata {
class symbolic_state_collection;
typedef boost::shared_ptr<symbolic_state_collection> symbolic_state_collection_ptr;
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
class reachability_scenario;
}

namespace hybrid_automata {

/** Parse the string s, instantiate a symbolic_state_collection in which the
 * sets of states are given by predicates and then adapt it according to the
 * scenario scen.
 *
 * The symbols used in s can be prefixed by a context. E.g., if the user
 * writes x,y,z to refer to variables in a component called sys, one
 * can provide the context "sys", so variables will be searched within
 * that component. Practically, the prefix "sys." is added to achieve this.
 * The variables are internally called sys.x,sys.y,sys.z. */
symbolic_state_collection_ptr parse_symbolic_state_collection_with_scenario(
		const reachability_scenario& scen, const std::string& s, const std::string& context="",
		const parser::parse_policy& ppol = parser::parse_policy());

/** Parse and adapt all automata in the file file_name.
 *
 * If sys_name is not empty, only this automaton is actually instantiated. */
std::vector<hybrid_automata::hybrid_automaton_ptr> parse_models_in_file_with_scenario(const reachability_scenario& scen,
		const std::string& file_name, bool adapt_all=true, const std::string& sys_name = "");

/** Parse the automata in file file_name and adapt and return the automaton
 * with name "system". */
hybrid_automaton_ptr parse_model_file_with_scenario(const reachability_scenario& scen,
		const std::string& file_name, const std::string& sys_name = "");

}

#endif /* PARSE_WITH_SCENARIO_H_ */
