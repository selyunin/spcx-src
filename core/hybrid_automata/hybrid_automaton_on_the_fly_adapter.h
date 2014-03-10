/*
 * hybrid_automaton_on_the_fly_adapter.h
 *
 *  Created on: Jun 24, 2013
 *      Author: frehse
 */

#ifndef HYBRID_AUTOMATON_ON_THE_FLY_ADAPTER_H_
#define HYBRID_AUTOMATON_ON_THE_FLY_ADAPTER_H_

#include "core/hybrid_automata/hybrid_automaton_wrapper.h"

namespace hybrid_automata {

/** Forward declarations */
class location;
typedef boost::shared_ptr<location> location_ptr;

/** This is an adapter class for locations
 *
 * It does not use a visitor pattern for more flexibility.
 */
class location_adapter {
public:
	typedef boost::shared_ptr<location_adapter> ptr;
	virtual ~location_adapter() {
	}
	;
	/** Default behavior is to do return the location without any modification. */
	virtual location_ptr accept(hybrid_automaton::ptr aut, location_ptr l) {
		return l;
	}
	;
};

/** This is an adapter class for transitions
 *
 * It does not use a visitor pattern for more flexibility.
 */
class transition_adapter {
public:
	typedef boost::shared_ptr<transition_adapter> ptr;
	virtual ~transition_adapter() {
	}
	;
	/** Default behavior is to do return the transition without any modification. */
	virtual transition_ptr accept(hybrid_automaton::ptr aut, transition_ptr t) {
		return t;
	}
	;
};

/** This class takes an hybrid automaton implementation and a set of adapters; whenever a location or
 *  a transition is requested, it is run through the adapters first.
 *
 *  @TODO: The result is cached in an explicit_automaton.
 *
 *  The adapters are passed as objects so that they can be given state information if necessary.
 *  Shared pointers are used to avoid problems with memory management -- this way there is no
 *  confusion about ownership etc.
 */
class hybrid_automaton_on_the_fly_adapter: public hybrid_automaton_wrapper {
public:
	/** Create an adapter */
	hybrid_automaton_on_the_fly_adapter(location_adapter::ptr lad,
			transition_adapter::ptr tad) {
		my_lad = lad;
		my_tad = tad;
	}

	/** Overloaded retrieval of locations */
	location_ptr get_location(const location_id& id) const {
		hybrid_automaton_on_the_fly_adapter* nonconst_this = const_cast<hybrid_automaton_on_the_fly_adapter*>(this);

		location_ptr orig_loc = hybrid_automaton_wrapper::get_location(id);
		if (my_lad) {
			return my_lad->accept(nonconst_this->get_ptr(),orig_loc);
		} else {
			return orig_loc;
		}
	}

	/** Overloaded retrieval of transitions */
	virtual transition_ptr get_transition(const transition_id& id) const {
		hybrid_automaton_on_the_fly_adapter* nonconst_this = const_cast<hybrid_automaton_on_the_fly_adapter*>(this);

		transition_ptr orig_trans = hybrid_automaton_wrapper::get_transition(
				id);
		if (my_tad) {
			return my_tad->accept(nonconst_this->get_ptr(),orig_trans);
		} else {
			return orig_trans;
		}
	}

private:
	location_adapter::ptr my_lad;
	transition_adapter::ptr my_tad;
};

}

#endif /* HYBRID_AUTOMATON_ON_THE_FLY_ADAPTER_H_ */
