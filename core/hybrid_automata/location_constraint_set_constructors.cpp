/*
 * location_constraint_set_constructors.cpp
 *
 *  Created on: Dec 29, 2009
 *      Author: frehse
 */

#include "location_constraint_set_constructors.h"
#include "location_eq_node.h"
#include "core/discrete/singleton_set.h"

namespace hybrid_automata {

class lcs_visitor: public tree::node_visitor {
public:

	lcs_visitor(hybrid_automata::location_constraint_set& lcs) :
		my_lcs(lcs) {
	}
	;

	/** By default, all nodes throw. We define only the
	 * ones we expect to see: AND nodes and location_eq:neq_nodes. */
	void visit(const valuation_functions::boolean_node* p) {
		if (p->my_op != AND)
			throw std::runtime_error("unacceptable boolean node in lcs_visitor");
		p->child1->accept(*this);
		p->child2->accept(*this);
	}
	;

	void visit(const hybrid_automata::location_eq_node* p) {
		if(p->get_equal()){
			hybrid_automata::location_constraint con(p->get_location_id(), true);
			my_lcs.add_constraint(p->get_automaton_id(), con);
		}
		else{
			hybrid_automata::location_constraint con(p->get_location_id(), false);
			my_lcs.add_constraint(p->get_automaton_id(), con);
		}
	}
	;

private:
	hybrid_automata::location_constraint_set& my_lcs;
};

hybrid_automata::location_constraint_set construct_location_constraint_set(
		const tree::node::ptr& p) {
	hybrid_automata::location_constraint_set lcs;
	 // if null, return the universe set
	if (p) {
		lcs_visitor v(lcs);
		p->accept(v);
	}
	return lcs;
}

discrete::discrete_set_ptr construct_discrete_set_from_conjunction(const tree::node::ptr& p) {
	// create the discrete set
	hybrid_automata::location_constraint_set lcs = construct_location_constraint_set(p);
	discrete::singleton_set::ptr d = discrete::singleton_set::ptr(new discrete::singleton_set(lcs));
	return d;
}

}

