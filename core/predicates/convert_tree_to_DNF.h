/*
 * convert_tree_to_DNF.h
 *
 *  Created on: Aug 26, 2009
 *      Author: gvincent
 */

#ifndef CONVERT_TREE_TO_DNF_H_
#define CONVERT_TREE_TO_DNF_H_

#include "utility/tree_traverser.h"
#include "core/predicates/valuation_function_tree_nodes.h"

namespace valuation_functions {

class has_OR_operator : public tree::tree_traverser {
public:
	has_OR_operator() {
		operator_or = false;
	}

	bool eval(tree::node::ptr p) {
		if (boolean_node* q = dynamic_cast<boolean_node*>(p.get())) {
			if (q->my_op == OR) {
				operator_or = true;
				return false;
			}
		}
		return true;
	}
	bool operator_or;
};

class DNF_test : public tree::tree_traverser {
public:
	DNF_test() {
		dnf = true;
	}

	bool eval(tree::node::ptr p) {
		if (boolean_node* q = dynamic_cast<boolean_node*>(p.get())) {
			if (q->my_op == AND) {
				has_OR_operator g;
				g.traverse_preorder(p);
				if(g.operator_or == true){
					dnf = false;
					return false;
				}
			}
			else if (q->my_op == NOT) {
				if (dynamic_cast<boolean_node*>(q->child1.get())) {
					dnf = false;
					return false;
				}
			}
		}
		return true;
	}
	bool dnf;
};

/**
 * Use class DNF_test to determine if the tree is a DNF or not
 */
bool is_DNF(const tree::node::ptr& p);

/**
 * check if p is boolean_node with operator OR
 */
bool is_OR_node(const tree::node::ptr& p);

/**
 * Create a node OR with children AND (to convert a node AND with children OR to DNF)
 * @param child1 boolean_node with operator OR
 */
tree::node::ptr create_OR_AND_AND(const boolean_node* child1, const tree::node::ptr& child2);

/**
 * if up_OR_operator find a AND node with child OR, up the OR to convert tree to DNF form
 */
tree::node::ptr convert_tree_to_DNF(const tree::node::ptr& p);

tree::node::ptr up_or_operator(const tree::node::ptr& p);

tree::node::ptr down_not_operator(const tree::node::ptr& p);


}

#endif /* CONVERT_TREE_TO_DNF_H_ */
