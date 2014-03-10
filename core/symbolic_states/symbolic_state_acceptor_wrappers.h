/*
 * symbolic_state_acceptor_wrappers.h
 *
 *  Created on: Sep 26, 2009
 *      Author: frehse
 */

#ifndef SYMBOLIC_STATE_ACCEPTOR_WRAPPERS_H_
#define SYMBOLIC_STATE_ACCEPTOR_WRAPPERS_H_

#include "boost/shared_ptr.hpp"
#include "core/symbolic_states/symbolic_state_acceptor.h"

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
typedef boost::shared_ptr<const hybrid_automaton> hybrid_automaton_const_ptr;
class symbolic_state;
typedef boost::shared_ptr<symbolic_state> symbolic_state_ptr;
class symbolic_state_collection;
typedef boost::shared_ptr<symbolic_state_collection> symbolic_state_collection_ptr;
typedef boost::shared_ptr<const symbolic_state_collection> symbolic_state_collection_const_ptr;
class post_operator;
typedef boost::shared_ptr<post_operator> post_operator_ptr;
class passed_and_waiting_list;
typedef boost::shared_ptr<passed_and_waiting_list> passed_and_waiting_list_ptr;
class reachability_scenario;
}

namespace hybrid_automata {

/** Wrapper classes to pass post operators as filters. */
class plwl_sink: public symbolic_state_acceptor {
public:
	plwl_sink(passed_and_waiting_list& plwl,hybrid_automaton_ptr& aut);
	/** Process the symbolic state sstate. */
	void accept(const symbolic_state_ptr& sstate);
private:
	passed_and_waiting_list& my_plwl;
	hybrid_automaton_ptr& my_aut;
};

class symbolic_state_collection_sink: public symbolic_state_acceptor {
public:
	symbolic_state_collection_sink(symbolic_state_collection& scoll);
	/** Process the symbolic state sstate. */
	void accept(const symbolic_state_ptr& sstate);
private:
	symbolic_state_collection& my_scoll;
};

}

#endif /* SYMBOLIC_STATE_ACCEPTOR_WRAPPERS_H_ */
