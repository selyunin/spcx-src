#ifndef HYBRID_AUTOMATON_NETWORK_H_
#define HYBRID_AUTOMATON_NETWORK_H_

#include "core/hybrid_automata/hybrid_automaton.h"

namespace hybrid_automata {

/** \brief An interface to handle networks of hybrid automata.
 *
 * This interface provides functions to add automata and a pointer to
 * an automaton that represents their parallel composition.
 */

class hybrid_automaton_network : public hybrid_automaton {

public:
	typedef boost::shared_ptr<hybrid_automaton_network> ptr;
	typedef boost::shared_ptr<const hybrid_automaton_network> const_ptr;

	virtual ~hybrid_automaton_network() {
	}
	;

	/** Virtual constructor for networks. */
	virtual hybrid_automaton_network* create_network() const = 0;

	/** \brief Return a pointer to a network composed with aut. The pointer
	 * can be *this. */
	virtual hybrid_automaton_network::ptr
			compute_or_assign_composition(hybrid_automaton::ptr aut) = 0;

	/** Obtain the ids of the automata in the network. */
	virtual automaton_id_set get_automata() const = 0;
};

}

#endif /*HYBRID_AUTOMATON_NETWORK_H_*/
