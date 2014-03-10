#ifndef AUTOMATON_SCENARIO_H_
#define AUTOMATON_SCENARIO_H_

#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/hybrid_automaton_network.h"

namespace hybrid_automata {

/** An automaton_scenario can construct an automaton and an automaton
 * network, and a representation of sets of states. */

class automaton_scenario {
public:
	virtual ~automaton_scenario() {};
	virtual hybrid_automata::hybrid_automaton::ptr create_automaton() = 0;
	virtual hybrid_automata::hybrid_automaton_network::ptr create_automaton_network() = 0;
};

}

#endif /*AUTOMATON_SCENARIO_H_*/
