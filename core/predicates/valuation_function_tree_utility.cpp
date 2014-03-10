/*
 * valuation_function_tree_utility.cpp
 *
 *  Created on: Aug 31, 2009
 *      Author: frehse
 */

#include "core/predicates/valuation_function_tree_utility.h"

#include "core/predicates/valuation_function_tree_nodes.h"
#include "utility/tree_traverser.h"

#include "io/common_input/scalar_node_creation.h" // for creating the node zero

namespace valuation_functions {

tree::node_ptr boolean_and(const tree::node_ptr& p1, const tree::node_ptr& p2) {
	if(p1) {
		if (p2) {
			return tree::node_ptr(new valuation_functions::boolean_node(AND, p1, p2));
		} else
			return p1;
	}
	else
		return p2;
}

class variable_id_getter: public tree::const_tree_traverser {
public:
	variable_id_getter();

	bool eval(tree::node::const_ptr p);
	variable_id_set vis;
};

variable_id_getter::variable_id_getter() {
}

bool variable_id_getter::eval(tree::node::const_ptr p) {
	if (const variable_node * q = dynamic_cast<const variable_node*> (p.get())) {
		vis.insert(q->my_id);
	}
	return true;
}

variable_id_set get_variable_ids(tree::node::const_ptr p) {
	variable_id_getter g;
	g.traverse_preorder(p);
	return g.vis;
}

class unprimed_variable_id_getter: public tree::tree_traverser {
public:
	unprimed_variable_id_getter();

	bool eval(tree::node_ptr p);
	variable_id_set vis;
};

unprimed_variable_id_getter::unprimed_variable_id_getter() {
}

bool unprimed_variable_id_getter::eval(tree::node::ptr p) {
	if (variable_node * q = dynamic_cast<variable_node*> (p.get())) {
		vis.insert(variable::get_primed_id(q->my_id, 0));
	}
	return true;
}

variable_id_set get_unprimed_variable_ids(tree::node::ptr p) {
	unprimed_variable_id_getter g;
	g.traverse_preorder(p);
	return g.vis;
}

class primed_variable_id_getter: public tree::tree_traverser {
public:
	primed_variable_id_getter(unsigned int pc) :
		my_prime_count(pc) {
	}
	;

	bool eval(tree::node::ptr p) {
		if (variable_node * q = dynamic_cast<variable_node*> (p.get())) {
			if (variable::get_prime_count(q->my_id) == my_prime_count) {
				vis.insert(variable::get_primed_id(q->my_id, 0));
			}
		}
		return true;
	}
	;
	variable_id_set vis;
private:
	unsigned int my_prime_count;
};

variable_id_set get_primed_variable_ids(tree::node_ptr p, unsigned int prime_count) {
	primed_variable_id_getter g(prime_count);
	g.traverse_preorder(p);
//std::cerr << g.vis;
	return g.vis;
}
;

class increase_primedness_traverser: public tree::tree_traverser {
public:
	increase_primedness_traverser(unsigned int delta) :
		my_delta(delta) {
	}
	;

	bool eval(tree::node_ptr p) {
		if (variable_node * q = dynamic_cast<variable_node*> (p.get())) {
			q->my_id = variable::get_primed_id(q->my_id, my_delta);
		}
		return true;
	}
private:
	unsigned int my_delta;
};

void increase_primedness(tree::node_ptr p, unsigned int delta) {
	increase_primedness_traverser t(delta);
	t.traverse_preorder(p);
}

void add_const_relation_constraints(tree::node_ptr& p,
		const variable_id_set& vis) {
	for (variable_id_set::const_iterator vit = vis.begin(); vit != vis.end(); ++vit) {
		variable_id v_id = *vit;

		variable_id primed_id = variable::get_id_primedness_increased(v_id);

		tree::node::ptr pred_left = tree::node::ptr(
				new valuation_functions::variable_node(primed_id));
		tree::node::ptr pred_right = tree::node::ptr(
				new valuation_functions::variable_node(v_id));
		tree::node::ptr pred_op = tree::node_ptr(
				new valuation_functions::comparison_node(EQ, pred_left,
						pred_right));

		p = boolean_and(p, pred_op);
	}

}

void add_const_dynamics_constraints(tree::node_ptr& p,
		const variable_id_set& vis) {
	for (variable_id_set::const_iterator vit = vis.begin(); vit != vis.end(); ++vit) {
		variable_id v_id = *vit;
		variable_id primed_id = variable::get_id_primedness_increased(v_id);

		tree::node::ptr pred_left = tree::node::ptr(
				new valuation_functions::variable_node(primed_id));
		tree::node::ptr pred_right =
				predicate_parser::create_scalar_const_node("0");
		tree::node::ptr pred_op = tree::node_ptr(
				new valuation_functions::comparison_node(EQ, pred_left,
						pred_right));

		p = boolean_and(p, pred_op);
	}
}

}

