#ifndef hybrid_automatON_VISITOR_H_
#define hybrid_automatON_VISITOR_H_

#include "core/hybrid_automata/location_id.h"
#include "core/hybrid_automata/transition_id.h"
#include "core/hybrid_automata/hybrid_automaton_network.h"
#include <iostream> /* for std::ostream -wsc */

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
class location;
class transition;
}

namespace hybrid_automata {

class hybrid_automaton_visitor {
public:
	typedef enum { locations_before_transitions, outgoing_transitions }  visiting_order;
	typedef enum { all_components_and_subsets, all_base_components, only_composition } network_visits;
	virtual ~hybrid_automaton_visitor();
	virtual void prologue(hybrid_automaton& h);
	virtual void visit(hybrid_automaton& h, location& l, location_id l_id);
	virtual void visit(hybrid_automaton& h, transition& l, transition_id t_id);
	virtual void epilogue(hybrid_automaton& e);
	virtual visiting_order get_visiting_order();
	/** Determines whether visitor should be passed on to all network components (for adapting dynamics, etc.) or
	 *  to composed automaton only (for printing).
	 *
	 * @note If only_composition, this means that the visitor is passed only to the composition generated so far
	 * when using an on-the-fly construction.
	 */
	virtual network_visits get_network_visits();
};

/** Output the automaton as stream. */
class print_visitor: public hybrid_automaton_visitor {
public:
	print_visitor(std::ostream& new_os);
	virtual ~print_visitor();
	virtual void prologue(hybrid_automaton& h);
	virtual void epilogue(hybrid_automaton& h);
	virtual void visit(hybrid_automaton& h, location& l, location_id l_id = 0);
	virtual void visit(hybrid_automaton& h, transition& t, transition_id t_id = 0);
	virtual network_visits get_network_visits();
private:
	std::ostream& os;
};

}

#endif /*hybrid_automatON_VISITOR_H_*/
