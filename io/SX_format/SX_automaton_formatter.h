/*
 * SX_automaton_formatter.h
 *
 *  Created on: Sep 16, 2010
 *      Author: frehse
 */

#ifndef SX_automaton_formatter_H_
#define SX_automaton_formatter_H_

#include "core/hybrid_automata/hybrid_automaton_visitor.h"

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
class location;
class transition;
}


namespace io {

/** Output the automaton as stream. */
class SX_automaton_formatter: public hybrid_automata::hybrid_automaton_visitor {
public:
	SX_automaton_formatter(std::ostream& new_os);
	virtual ~SX_automaton_formatter();
	virtual void file_prologue();
	virtual void file_epilogue();
	virtual std::string format_ident(const std::string& ident,const hybrid_automata::hybrid_automaton& h);
	virtual void prologue(hybrid_automata::hybrid_automaton& h);
	virtual void epilogue(hybrid_automata::hybrid_automaton& h);
	virtual void visit(hybrid_automata::hybrid_automaton& h, hybrid_automata::location& l, hybrid_automata::location_id l_id = 0);
	virtual void visit(hybrid_automata::hybrid_automaton& h, hybrid_automata::transition& t, hybrid_automata::transition_id t_id = 0);
	virtual network_visits get_network_visits();
private:
	std::ostream& os;
};

}


#endif /* SX_automaton_formatter_H_ */
