/*
 * node_eval_visitor.h
 *
 *  Created on: Aug 28, 2009
 *      Author: gvincent
 */

#ifndef NODE_EVAL_VISITOR_H_
#define NODE_EVAL_VISITOR_H_

#include <math.h>
#include "core/predicates/valuation_function_tree_nodes.h"
#include "core/predicates/valuation.h"
#include "core/hybrid_automata/location_tree.h"
#include "math/type_conversion.h"
#include "math/vdom/vdom_vector.h"
//#include "../utility/shared_ptr_output.h"


namespace valuation_functions {

using ::sqrt;
using std::sqrt;
using ::pow;
using std::pow;

template<typename bool_type, typename scalar_type>
class node_eval_visitor: public tree::node_visitor {
public:

	/**
	 * Four constructor:
	 * first with variable_valuation (only arithmetic expressions)
	 * second with automaton_valuation (only state expressions)
	 * third with both
	 * the last has no valuation
	 */
	node_eval_visitor(const math::vdom_vector<scalar_type>& val) :
		my_valuation(val) {
		my_bool = bool_type(true);
		my_scalar = scalar_type(0);
		my_automaton_valuation
				= hybrid_automata::automaton_valuation<hybrid_automata::location_id>::ptr(
						new hybrid_automata::automaton_valuation<hybrid_automata::location_id>());
	}
	;

	node_eval_visitor(const typename hybrid_automata::automaton_valuation<
			hybrid_automata::location_id>::const_ptr& automaton_val) :
		my_automaton_valuation(automaton_val) {
		my_bool = bool_type(true);
		my_scalar = scalar_type(0);
		my_valuation = math::vdom_vector<scalar_type>::const_reference(new math::vdom_vector<scalar_type>());
	}
	;

	node_eval_visitor(
			const math::vdom_vector<scalar_type>& val,
			const typename hybrid_automata::automaton_valuation<hybrid_automata::location_id>::const_ptr& automaton_val) :
		my_valuation(val), my_automaton_valuation(automaton_val) {
		my_bool = bool_type(true);
		my_scalar = scalar_type(0);
	}
	;

	node_eval_visitor() {
		my_bool = bool_type(true);
		my_scalar = scalar_type(0);
		my_valuation = math::vdom_vector<scalar_type>::const_reference(new math::vdom_vector<scalar_type>());
		my_automaton_valuation
				= hybrid_automata::automaton_valuation<hybrid_automata::location_id>::ptr(
						new hybrid_automata::automaton_valuation<hybrid_automata::location_id>());
	}
	;

	virtual ~node_eval_visitor() {
	}
	;

	/**
	 * Like all node_visitor, node_eval_visitor accepts five const type:
	 * Rational, bool, calc_string, int and double.
	 * Convert each type to scalar_type with convert_element.
	 */
	virtual void visit(const const_node<Rational>* p) {
		assert(p);
		my_scalar = convert_element<scalar_type, Rational> (p->my_val);
	}
	;
	virtual void visit(const const_node<bool>* p) {
		assert(p);
		my_scalar = convert_element<scalar_type, bool> (p->my_val);
	}
	;
	virtual void visit(const const_node<calc_string>* p) {
		assert(p);
		my_scalar = convert_element<scalar_type, calc_string> (p->my_val);
	}
	;
	virtual void visit(const const_node<int>* p) {
		assert(p);
		my_scalar = convert_element<scalar_type, int> (p->my_val);
	}
	;
	virtual void visit(const const_node<double>* p) {
		assert(p);
		my_scalar = convert_element<scalar_type, double> (p->my_val);
	}
	;

	virtual void visit(const variable_node* p) {
		assert(p);
		my_scalar = my_valuation.get_coeff_with_id(p->my_id);
	}
	;

	/**
	 * Evaluate first child, and store its arithmetic value in scalar_type temp.
	 * Then, evaluate second child and perform the operation.
	 * Store the result in my_scalar.
	 */
	virtual void visit(const arithmetic_node* p) {
		assert(p);
		if(p->my_op == NEG){
			p->child1->accept(*this);
			my_scalar = -my_scalar;
		}
		else if(p->my_op == SQRT)
		{
			p->child1->accept(*this);
			my_scalar = sqrt(my_scalar);
		}
		else if(p->my_op == LOG)
		{
			p->child1->accept(*this);
			my_scalar = log(my_scalar); // natural logarithm
		}
		else if(p->my_op == SIN)
		{
			p->child1->accept(*this);
			my_scalar = sin(my_scalar); // input angle in rad
		}
		else if(p->my_op == EXP)
		{
			p->child1->accept(*this);
			my_scalar = exp(my_scalar);
		}
		else if(p->my_op == TAN)
		{
			p->child1->accept(*this);
			my_scalar = tan(my_scalar);
		}
		else if(p->my_op == COS)
		{
			p->child1->accept(*this);
			my_scalar = cos(my_scalar);
		}
		else {
			p->child1->accept(*this);
			scalar_type temp = my_scalar;
			p->child2->accept(*this);

			switch (p->my_op) {
			case ADD:
				my_scalar = temp + my_scalar;
				break;
			case SUB:
				my_scalar = temp - my_scalar;
				break;
			case MUL:
				my_scalar = temp * my_scalar;
				break;
			case DIV:
				my_scalar = temp / my_scalar;
				break;
			case POW:
				my_scalar = pow(temp, my_scalar);
				break;
			default:
				throw std::runtime_error("unknown arithmetic operation");
				break;
			}
		}
	}
	;

	/**
	 * Evaluate first child, and store its boolean value in bool_type temp.
	 * Then, evaluate second child and apply the operator.
	 * Store the result in my_bool.
	 */
	virtual void visit(const boolean_node* p) {
		assert(p);
		if (p->my_op != NOT) {
			p->child1->accept(*this);
			bool_type temp = my_bool;
			p->child2->accept(*this);

			if (p->my_op == AND) {
				my_bool = (temp && my_bool);
			} else if (p->my_op == OR) {
				my_bool = (temp || my_bool);
			} else
				throw std::runtime_error("unknown boolean operation");

		} else {
			p->child1->accept(*this);
			my_bool = !my_bool;
		}
	}
	;

	/**
	 * Evaluate first child, and store its arithmetic value in scalar_type temp.
	 * Then, evaluate second child and compare both.
	 * Store the result in my_bool.
	 */
	virtual void visit(const comparison_node* p) {
		assert(p);
		p->child1->accept(*this);
		scalar_type temp = my_scalar;
		p->child2->accept(*this);

		switch (p->my_op) {
		case EQ:
			my_bool = (temp == my_scalar);
			break;
		case LE:
			my_bool = (temp <= my_scalar);
			break;
		case LT:
			my_bool = (temp < my_scalar);
			break;
		case GE:
			my_bool = (temp >= my_scalar);
			break;
		case GT:
			my_bool = (temp > my_scalar);
			break;
		default:
			throw std::runtime_error("unknown comparison operation");
			break;
		}

	}
	;

	/**
	 * Evaluate state expression with my_automaton_valuation.
	 * Store the result in my_bool.
	 */
	virtual void visit(const hybrid_automata::location_eq_node* p) {
		assert(p);
		if (p->get_equal())
			my_bool = (bool_type) (p->get_location_id() == my_automaton_valuation->get_value(
					p->get_automaton_id()));
		else
			my_bool = (bool_type) (p->get_location_id() != my_automaton_valuation->get_value(
					p->get_automaton_id()));
	}
	;

	bool_type get_bool() {
		return my_bool;
	}

	scalar_type get_scalar() {
		return my_scalar;
	}

private:
	bool_type my_bool;
	scalar_type my_scalar;

	typename math::vdom_vector<scalar_type> my_valuation;
	typename hybrid_automata::automaton_valuation<hybrid_automata::location_id>::const_ptr
			my_automaton_valuation;
};

/** Evaluates an arithmetic expression given by a tree.
 * The expressions are evaluated as objects of type scalar_type.
 * Variable values are stored in a variable_valuation.
 */
template<typename scalar_type>
scalar_type arithmetic_eval_node(const tree::node::ptr& p, const typename math::vdom_vector<scalar_type>& val) {

	if (p) {
		node_eval_visitor<bool, scalar_type> v(val);
		p->accept(v);
		return v.get_scalar();

	} else
		return (scalar_type) (0);
}

/** Evaluates an boolean expression given by a tree.
 * Variable values are stored in a variable_valuation.
 */
template<typename bool_type, typename scalar_type>
bool_type boolean_eval_node(const tree::node::ptr& p,
		const typename math::vdom_vector<scalar_type>& val) {
	if (p) {
		node_eval_visitor<bool_type, scalar_type> v(val);
		p->accept(v);
		return v.get_bool();
	} else
		return (bool_type) (true);
}

/** Evaluates an state expression (boolean expression combining arithmetic and state expressions) given by a tree.
 * Variable values are stored in a variable_valuation and states are stored in a automaton_valuation.
 */
template<typename bool_type, typename scalar_type>
bool_type state_eval_node(
		const tree::node::ptr& p,
		const typename math::vdom_vector<scalar_type>& val,
		const typename hybrid_automata::automaton_valuation<hybrid_automata::location_id>::const_ptr& automaton_val) {
	if (p) {
		node_eval_visitor<bool_type, scalar_type> v(val, automaton_val);
		p->accept(v);
		return v.get_bool();
	} else
		return (bool_type) (true);
}

}

#endif /* NODE_EVAL_VISITOR_H_ */
