/*
 * tree_node_visitor.h
 *
 *  Created on: Aug 28, 2009
 *      Author: gvincent
 */

#ifndef TREE_NODE_VISITOR_H_
#define TREE_NODE_VISITOR_H_

#include <stdexcept>
#include "utility/calc_string.h"
#include "math/scalar_types/rational.h"
#include "math/matrix.h"

namespace tree {

class node_visitor {
public:

	node_visitor() {
	}
	;

	virtual ~node_visitor() {
	}
	;

	virtual void visit(const valuation_functions::const_node<Rational>* p) {
		throw std::runtime_error("Missing implementation");
	}
	;

	virtual void visit(const valuation_functions::const_node<math::matrix<Rational> >* p) {
			throw std::runtime_error("Missing implementation");
		}
		;

	virtual void visit(const valuation_functions::const_node<bool>* p) {
		throw std::runtime_error("Missing implementation");
	}
	;

	virtual void visit(const valuation_functions::const_node<calc_string>* p) {
		throw std::runtime_error("Missing implementation");
	}
	;

	virtual void visit(const valuation_functions::const_node<int>* p) {
		throw std::runtime_error("Missing implementation");
	}
	;

	virtual void visit(const valuation_functions::const_node<double>* p) {
		throw std::runtime_error("Missing implementation");
	}
	;

	virtual void visit(const valuation_functions::const_node<long double>* p) {
		throw std::runtime_error("Missing implementation");
	}
	;

	virtual void visit(const valuation_functions::variable_node* p) {
		throw std::runtime_error("Missing implementation");
	}
	;
	virtual void visit(const valuation_functions::arithmetic_node* p) {
		throw std::runtime_error("Missing implementation");
	}
	;
	virtual void visit(const valuation_functions::boolean_node* p) {
		throw std::runtime_error("Missing implementation");
	}
	;
	virtual void visit(const valuation_functions::comparison_node* p) {
		throw std::runtime_error("Missing implementation");
	}
	;
	virtual void visit(const hybrid_automata::location_eq_node* p) {
		throw std::runtime_error("Missing implementation");
	}
	;
};

}

#endif /* TREE_NODE_VISITOR_H_ */
