/*
 * node_print_visitor.cpp
 *
 *  Created on: Sep 2, 2009
 *      Author: frehse
 */

#include "core/predicates/node_print_visitor.h"

#include <sstream>
//#include "../parser/state_parser.h"
#include "core/hybrid_automata/automaton_cache.h"
#include "core/hybrid_automata/hybrid_automaton.h"
#include "core/hybrid_automata/location.h"
#include "core/hybrid_automata/location_eq_node.h"
#include "utility/calc_string_operators.h"

namespace valuation_functions {

/** Returns for a binary operator whether there should be parentheses around
 * the first and the second operand.
 */
std::pair<bool, bool> needs_parentheses(const arithmetic_node* p) {
	assert(p);
	bool par1 = false;
	bool par2 = false;
	if (p->my_op == MUL || p->my_op == DIV) {
		if (arithmetic_node * q = dynamic_cast<arithmetic_node*> (p->child1.get())) {
			if (q->my_op == ADD || q->my_op == SUB) {
				par1 = true;
			}
		}
		if (arithmetic_node * q = dynamic_cast<arithmetic_node*> (p->child2.get())) {
			if (q->my_op == ADD || q->my_op == SUB || q->my_op == NEG) {
				par2 = true;
			}
		}
	} else if (p->my_op == NEG || p->my_op == POW) {
		if (dynamic_cast<tree::binary_node*> (p->child1.get())) {
			par1 = true;
		}
	}
	else if (p->my_op == SQRT || p->my_op == LOG || p->my_op == SIN || p->my_op == EXP
			|| p->my_op == TAN || p->my_op == COS) {

			par1 = true;

	}
	return std::make_pair(par1, par2);
//	return std::make_pair(true, true);
}

/** Returns for a binary operator whether there should be parentheses around
 * the first and the second operand.
 */
std::pair<bool, bool> needs_parentheses(const boolean_node* p) {
	assert(p);
	bool par1 = false;
	bool par2 = false;
	if (p->my_op == AND) {
		if (boolean_node * q = dynamic_cast<boolean_node*> (p->child1.get())) {
			if (q->my_op == OR) {
				par1 = true;
			}
		}
		if (boolean_node * q = dynamic_cast<boolean_node*> (p->child2.get())) {
			if (q->my_op == OR) {
				par2 = true;
			}
		}
	} else if (p->my_op == NOT) {
		if (dynamic_cast<tree::binary_node*> (p->child1.get())) {
			par1 = true;
		}
	}
	return std::make_pair(par1, par2);
}

/** Returns for a binary operator whether there should be parentheses around
 * the first and the second operand.
 */
std::pair<bool, bool> needs_parentheses(const comparison_node* p) {
	assert(p);
	bool par1 = false;
	bool par2 = false;
	return std::make_pair(par1, par2);
}

std::string operator_to_string(const arithmetic_node* p) {
	if (p->my_op == NEG) {
		return "-";
	} else if (p->my_op == ADD) {
		return " + ";
	} else if (p->my_op == SUB) {
		return " - ";
	} else if (p->my_op == MUL) {
		return "*";
	} else if (p->my_op == DIV) {
		return "/";
	}
	else if (p->my_op == SQRT) {
		return "sqrt";
	}
	else if (p->my_op == SIN) {
		return "sin";
	}
	else if (p->my_op == LOG) {
		return "log";
	}
	else if (p->my_op == EXP) {
			return "exp";
	}
	else if (p->my_op == POW) {
		return "^";
	}
	else if (p->my_op == COS) {
		return "cos";
	}
	else if (p->my_op == TAN) {
		return "tan";
	}
	else {
		throw std::runtime_error("unknown arithmetic operation");
	}
}

std::string operator_to_string(const boolean_node* p) {
	if (p->my_op == NOT) {
		return "!";
	} else if (p->my_op == AND) {
		return " & ";
	} else if (p->my_op == OR) {
		return " | ";
	} else {
		throw std::runtime_error("unknown boolean operation");
	}
}

std::string operator_to_string(const comparison_node* p) {
	if (p->my_op == LT) {
		return " < ";
	} else if (p->my_op == LE) {
		return " <= ";
	} else if (p->my_op == GT) {
		return " > ";
	} else if (p->my_op == GE) {
		return " >= ";
	} else if (p->my_op == EQ) {
		return " == ";
	} else {
		throw std::runtime_error("unknown comparison operation");
	}
}

class node_print_visitor: public tree::node_visitor {
public:

	node_print_visitor(std::ostream& os) :
		my_os(os) {
	}

	virtual ~node_print_visitor(){
	}

	virtual void visit(const const_node<Rational>* p){
		assert(p);
		my_os << p->my_val;
	}
	virtual void visit(const const_node<math::matrix<Rational> >* p){
		assert(p);
		my_os << p->my_val;
	}
	virtual void visit(const const_node<bool>* p){
		assert(p);
		my_os << p->my_val;
	}
	virtual void visit(const const_node<calc_string>* p){
		assert(p);
		my_os << p->my_val;
	}
	virtual void visit(const const_node<int>* p){
		assert(p);
		my_os << p->my_val;
	}
	virtual void visit(const const_node<double>* p){
		assert(p);
		my_os << p->my_val;
	}

	virtual void visit(const variable_node* p){
		assert(p);
		my_os << variable(p->my_id);
	}

	virtual void visit(const arithmetic_node* p) {
		assert(p);
		std::pair<bool, bool> pars = needs_parentheses(p);
		if (p->my_op == NEG || p->my_op == SQRT || p->my_op == LOG || p->my_op == SIN || p->my_op == EXP
				|| p->my_op == TAN || p->my_op == COS) {
			my_os << operator_to_string(p);
				if (pars.first)
					my_os << "(";
				p->child1->accept(*this);
				if (pars.first)
					my_os << ")";

		} else {
			if (pars.first)
				my_os << "(";
			p->child1->accept(*this);
			if (pars.first)
				my_os << ")";
			my_os << operator_to_string(p);
			if (pars.second)
				my_os << "(";
			p->child2->accept(*this);
			if (pars.second)
				my_os << ")";

		}
		//my_os << "(" + v1.get_string() + arithmetic_operator_to_string(p->my_op) + v2.get_string() + ")";
	}

	virtual void visit(const boolean_node* p){
		assert(p);
		std::pair<bool, bool> pars = needs_parentheses(p);
		if (p->my_op != NOT) {
			if (pars.first)
				my_os << "(";
			p->child1->accept(*this);
			if (pars.first)
				my_os << ")";
			my_os << operator_to_string(p);
			if (pars.second)
				my_os << "(";
			p->child2->accept(*this);
			if (pars.second)
				my_os << ")";
		} else {
			my_os << operator_to_string(p);
			if (pars.first)
				my_os << "(";
			p->child1->accept(*this);
			if (pars.first)
				my_os << ")";
		}
		//my_os << "(" + v1.get_string() + boolean_operator_to_string(p->my_op) + v2.get_string() + ")";
	}

	virtual void visit(const comparison_node* p) {
		assert(p);
		std::pair<bool, bool> pars = needs_parentheses(p);
		if (pars.first)
			my_os << "(";
		p->child1->accept(*this);
		if (pars.first)
			my_os << ")";
		my_os << operator_to_string(p);
		if (pars.second)
			my_os << "(";
		p->child2->accept(*this);
		if (pars.second)
			my_os << ")";
		//my_os << "(" + v1.get_string() + comparison_operator_to_string(p->my_op) + v2.get_string() + ")";
	}

	virtual void visit(const hybrid_automata::location_eq_node* p) {
		assert(p);
		hybrid_automata::hybrid_automaton::ptr a =
				hybrid_automata::hybrid_automaton_cache::get_automaton(p->get_automaton_id());
		std::string aut = a->get_name();
		std::string loc;
		loc = a->get_location(p->get_location_id())->get_name();
		if(p->get_equal())
			my_os << "loc(" + aut + ") == " + loc + "";
		else
			my_os << "loc(" + aut + ") != " + loc + "";
	}

//	const std::string& get_string();

private:
	std::ostream& my_os;
};


void print_node(std::ostream& os, const tree::node* p) {
	if (p) {
		node_print_visitor v(os);
		p->accept(v);
	}
}

}

std::ostream& operator<<(std::ostream& os, const tree::node& p) {
	valuation_functions::node_print_visitor v(os);
	p.accept(v);
	return os;
}

std::ostream& operator<<(std::ostream& os, const tree::node::ptr& p) {
	if(p)
		valuation_functions::print_node(os, p.get());
	else
			os << "";
	return os;
}

std::ostream& print_as_predicate(std::ostream& os, const tree::node::ptr& p) {
	if (p)
		valuation_functions::print_node(os, p.get());
	else
		os << "true";
	return os;
}

std::ostream& print_as_arithmetic(std::ostream& os, const tree::node::ptr& p) {
	if (p)
		valuation_functions::print_node(os, p.get());
	else
		os << "0";
	return os;
}
