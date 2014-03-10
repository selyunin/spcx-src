/*
 * network_creation.h
 *
 *  Created on: Sep 7, 2009
 *      Author: frehse
 */

#ifndef NETWORK_CREATION_H_
#define NETWORK_CREATION_H_

#include "boost/shared_ptr.hpp"

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
typedef boost::shared_ptr<const hybrid_automaton> hybrid_automaton_const_ptr;
class hybrid_automaton_network;
typedef boost::shared_ptr<hybrid_automaton_network> hybrid_automaton_network_ptr;
typedef boost::shared_ptr<const hybrid_automaton_network> hybrid_automaton_network_const_ptr;
class passed_and_waiting_list;
typedef boost::shared_ptr<passed_and_waiting_list> passed_and_waiting_list_ptr;
typedef boost::shared_ptr<const passed_and_waiting_list> passed_and_waiting_list_const_ptr;
class symbolic_state;
typedef boost::shared_ptr<symbolic_state> symbolic_state_ptr;
typedef boost::shared_ptr<const symbolic_state> symbolic_state_const_ptr;
class symbolic_state_collection;
typedef boost::shared_ptr<symbolic_state_collection> symbolic_state_collection_ptr;
typedef boost::shared_ptr<const symbolic_state_collection> symbolic_state_collection_const_ptr;
class post_operator;
typedef boost::shared_ptr<post_operator> post_operator_ptr;
typedef boost::shared_ptr<const post_operator> post_operator_const_ptr;
}

namespace parser {
namespace automaton_parser {

hybrid_automata::hybrid_automaton_network_ptr create_automaton_network(const std::string& name);

hybrid_automata::hybrid_automaton_ptr create_automaton(const std::string& name);
}
}

#endif /* NETWORK_CREATION_H_ */
