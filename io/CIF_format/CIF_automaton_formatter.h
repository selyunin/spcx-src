/*
 * CIF_automaton_formatter.h
 *
 *  Created on: May 11, 2011
 *      Author: Goyal
 */

#ifndef CIF_AUTOMATON_FORMATTER_H_
#define CIF_AUTOMATON_FORMATTER_H_

#include "core/hybrid_automata/hybrid_automaton_visitor.h"

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
class location;
class transition;
}


namespace io {

/** Output the automaton as stream. */
class CIF_automaton_formatter: public hybrid_automata::hybrid_automaton_visitor {
public:
	typedef enum { inv_predicate, reset_map_predicate, guard_predicate }  predicate_type;
	CIF_automaton_formatter(std::ostream& new_os);
	virtual ~CIF_automaton_formatter();
	virtual void file_prologue();
	virtual void file_epilogue();
	virtual std::string format_ident(const std::string& ident,const hybrid_automata::hybrid_automaton& h);
	virtual void prologue(hybrid_automata::hybrid_automaton& h);
	virtual void epilogue(hybrid_automata::hybrid_automaton& h);
	virtual hybrid_automata::hybrid_automaton_visitor::visiting_order get_visiting_order();
	virtual void visit(hybrid_automata::hybrid_automaton& h, hybrid_automata::location& l, hybrid_automata::location_id l_id = 0);
	virtual void visit(hybrid_automata::hybrid_automaton& h, hybrid_automata::transition& t, hybrid_automata::transition_id t_id = 0);

	virtual std::string spaceex_to_cif(const std::string& str, predicate_type pred);
	virtual void output(hybrid_automata::hybrid_automaton& h);
	std::string format_with_context(const std::string& ident, const std::string& context);

private:
	std::ostream& os;
	};

}


#endif /* CIF_AUTOMATON_FORMATTER_H_ */
