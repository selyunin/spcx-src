#ifndef SHARE_TREE_H_
#define SHARE_TREE_H_

#include "math/vdom/variable.h"
#include "utility/stl_helper_functions.h"
#include "utility/tree_traverser.h"
#include "core/predicates/valuation_function_tree_nodes.h"

#include "core/predicates/valuation_function_tree_utility.h"
#include "core/predicates/node_print_visitor.h"

namespace valuation_functions {

/**
 * Use class has_primed to determine if the tree contains a primed variable
 */
bool has_primed_variable(const tree::node::ptr& p);

/**
 * Check the tree contains only AND as logic operator (not OR)
 */
bool is_only_and(const tree::node::ptr& p);

/**
 * make a predicate (assignment free) tree from a tree with assignment and predicate.
 * extract predicate from the first tree and join them
 */
template<typename bool_type> tree::node::ptr make_assignment_free_tree(const tree::node::ptr& p) {
	tree::node::ptr new_node;

	if (boolean_node* q = dynamic_cast<boolean_node*>(p.get())) {
		tree::node::ptr child1 = make_assignment_free_tree<bool_type> (q->child1);

		if (q->my_op == NOT) {
			if (child1 == tree::node::null_node())
				new_node = tree::node::ptr(tree::node::null_node());
			else
				new_node = tree::node::ptr(q);
		} else if (q->my_op == AND) {
			tree::node::ptr child2 = make_assignment_free_tree<bool_type> (q->child2);

			if (child1 == tree::node::null_node() && child2 == tree::node::null_node())
				new_node = tree::node::ptr(tree::node::null_node());

			else if (child1 == tree::node::null_node())
				new_node = child2;

			else if (child2 == tree::node::null_node())
				new_node = child1;

			else
				new_node = tree::node::ptr(new boolean_node(AND, child1, child2));
		}

	}
	else if (const_node<bool_type> * q = dynamic_cast<const_node<bool_type>*> (p.get())) {
		new_node = tree::node::ptr(new const_node<bool_type> (q->my_val));
	} else if (dynamic_cast<comparison_node*> (p.get())) {
		if (has_primed_variable(p)) {
			//it's a flow predicate
			new_node = tree::node::ptr(tree::node::null_node());
		} else {
			new_node = p;
		}
	} else if (p == tree::node::null_node()) {
		new_node = p;
	} else
		throw std::runtime_error("unmatched node types0");

	return new_node;
}

/**
 * make a assignment tree from a tree with assignment and predicate.
 * extract assignments (and flow predicate) from the first tree and join them
 */
template<typename bool_type> tree::node::ptr make_assignment_only_tree(const tree::node::ptr& p) {
	tree::node::ptr new_node;

	if (boolean_node* q = dynamic_cast<boolean_node*>(p.get())) {
		tree::node::ptr child1 = make_assignment_only_tree<bool_type> (q->child1);

		if (q->my_op == NOT) {
			if (child1 == tree::node::null_node())
				new_node = tree::node::ptr(tree::node::null_node());
			else
				new_node = tree::node::ptr(q);
		} else if (q->my_op == AND) {
			tree::node::ptr child2 = make_assignment_only_tree<bool_type> (q->child2);

			if (child1 == tree::node::null_node() && child2 == tree::node::null_node())
				new_node = tree::node::ptr(tree::node::null_node());

			else if (child1 == tree::node::null_node())
				new_node = child2;

			else if (child2 == tree::node::null_node())
				new_node = child1;

			else
				new_node = tree::node::ptr(new boolean_node(AND, child1, child2));
		}

	}
	else if (dynamic_cast<const_node<bool_type>*> (p.get())) {
		new_node = tree::node::ptr(tree::node::null_node());
	} else if (dynamic_cast<comparison_node*> (p.get())) {
		if (!has_primed_variable(p))
			new_node = tree::node::ptr(tree::node::null_node());
		else
			//it's a flow predicate
			new_node = p;
	} else if (p == tree::node::null_node()) {
		new_node = p;
	} else
		throw std::runtime_error("unmatched node types1");

	return new_node;
}

/**
 * share a tree to two tree: a assignment free tree and a assignment only tree
 */
template<typename bool_type> std::pair<tree::node::ptr, tree::node::ptr> share_tree(
		tree::node::ptr p) {
	if (is_only_and(p)) {
		tree::node::ptr p1 = make_assignment_free_tree<bool_type> (p);
		tree::node::ptr p2 = make_assignment_only_tree<bool_type> (p);

		std::pair<tree::node::ptr, tree::node::ptr> two_tree = std::make_pair(p1, p2);
		return two_tree;
	} else
		throw std::runtime_error("Tree contains operator OR.");
}

/** Create a const relation over a set of variables
 *
 * Returns a predicate that is the conjunction of x'==x for all
 * x in vis.
 */
tree::node::ptr create_const_relation(const variable_id_set& vis);

}
#endif /*SHARE_TREE_H_*/
