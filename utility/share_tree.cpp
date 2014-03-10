/*
 * share_tree.cpp
 *
 *  Created on: Oct 17, 2010
 *      Author: gvincent
 */


#include "utility/share_tree.h"

namespace valuation_functions {

class has_primed: public tree::tree_traverser {
public:
	has_primed() {
		prime = false;
	}

	bool eval(tree::node::ptr p) {
		if (variable_node* q = dynamic_cast<variable_node*>(p.get())) {
			if (variable::get_prime_count(q->my_id) > 0) {
				prime = true;
				return false;
			}
		}
		return true;
	}
	bool prime;
};


bool has_primed_variable(const tree::node::ptr& p) {
	has_primed g;
	g.traverse_preorder(p);
	return g.prime;
}

bool is_only_and(const tree::node::ptr& p) {
	if (boolean_node* q = dynamic_cast<boolean_node*>(p.get())) {
		if (q->my_op == OR)
			return false;
		else if (q->my_op == NOT)
			return is_only_and(q->child1);
		else
			return is_only_and(q->child1) || is_only_and(q->child2);
	} else
		return true;
}

tree::node::ptr create_const_relation_constraint(const variable_id & vid) {
	tree::node::ptr var = tree::node::ptr(new variable_node(vid));
	tree::node::ptr pvar = tree::node::ptr(new variable_node(variable::get_primed_id(vid)));
	tree::node::ptr new_node = tree::node::ptr(new comparison_node(EQ, pvar, var));
	return new_node;
}
;

tree::node::ptr create_const_relation(const variable_id_set& vis) {

	tree::node::ptr root_node = tree::node::ptr();
	for (variable_id_set::const_iterator it = vis.begin(); it != vis.end(); ++it) {
		tree::node::ptr new_node = create_const_relation_constraint(*it);
		if (it == vis.begin()) {
			root_node = new_node;
		} else {
			root_node = tree::node::ptr(new boolean_node(AND, root_node, new_node));
		}
	}

	return root_node;
}

}
