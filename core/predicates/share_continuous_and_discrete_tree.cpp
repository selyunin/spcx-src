/*
 * share_continuous_and_discrete_tree.cpp
 *
 *  Created on: Sep 27, 2009
 *      Author: gvincent
 */

#include "io/common_input/state_parser.h"
#include "core/hybrid_automata/location_eq_node.h"
#include "core/predicates/convert_tree_to_DNF.h"
#include "core/predicates/share_continuous_and_discrete_tree.h"
#include "io/common_input/bool_node_creation.h"

namespace hybrid_automata {

inline tree::node::ptr equal_or_void(const tree::node::ptr& p, bool is_p) {
	if (is_p)
		return p;
	else
		return tree::node::null_node();
}

/**
 * If continuous == true, return a continuous tree (with all continuous nodes in p).
 * Else, return a discrete tree (with all discrete nodes in p).
 */
tree::node::ptr make_continuous_or_discrete_tree(
		const tree::node::ptr& p, bool continuous) {
	tree::node::ptr new_node;

	if (valuation_functions::boolean_node* q = dynamic_cast<valuation_functions::boolean_node*>(p.get())) {
		/**
		 * if node is AND, call itself for each child.
		 */
		if (q->my_op == AND) {
			tree::node::ptr child2 = make_continuous_or_discrete_tree(q->child2,
					continuous);
			tree::node::ptr child1 = make_continuous_or_discrete_tree(q->child1,
					continuous);

			/**
			 * if all children are void, node is void.
			 */
			if (child1 == tree::node::null_node() && child2 == tree::node::null_node())
				new_node = tree::node::null_node();

			/**
			 * if a child is void, node becomes the other child.
			 */
			else if (child1 == tree::node::null_node())
				new_node = child2;

			else if (child2 == tree::node::null_node())
				new_node = child1;

			else
				new_node = tree::node::ptr(new valuation_functions::boolean_node(AND, child1,
						child2));
		} else
			throw std::runtime_error("boolean node is not AND");

	} else if (dynamic_cast<location_eq_node*> (p.get()))
		new_node = equal_or_void(p, !continuous);

	else if (predicate_parser::is_bool_const_node(p))
		new_node = equal_or_void(p, continuous);

	else if (dynamic_cast<valuation_functions::comparison_node*> (p.get()))
		new_node = equal_or_void(p, continuous);

	else if (p == tree::node::null_node())
		new_node = p;
	else
		throw std::runtime_error("unmatched node types");

	return new_node;
}

/**
 * Share continuous and discrete trees and store them in a vector of pairs.
 * For all OR node, call itself for each child.
 * Else, make a continuous tree and a discrete tree with make_continuous_or_discrete_tree,
 * then compose a pair with them and store it in the vector.
 */
void share_continuous_and_discrete_tree_recursion(
		const tree::node::ptr& p,
		std::vector<std::pair<tree::node::ptr, tree::node::ptr> >& list_tree) {
	bool is_disjunction = false;
	if (valuation_functions::boolean_node* b = dynamic_cast<valuation_functions::boolean_node*>(p.get())) {
		if (b->my_op == OR) {
			share_continuous_and_discrete_tree_recursion(b->child1, list_tree);
			share_continuous_and_discrete_tree_recursion(b->child2, list_tree);
			is_disjunction = true;
		}
	}
	if (!is_disjunction) {
		list_tree.push_back(std::make_pair<tree::node::ptr, tree::node::ptr>(
				make_continuous_or_discrete_tree(p, true),
				make_continuous_or_discrete_tree(p, false)));
	}
}

/**
 * Convert p to DNF form, then use share_continuous_and_discrete_tree_recursion.
 */
void share_continuous_and_discrete_tree(const tree::node::ptr& p,
		std::vector<std::pair<tree::node::ptr, tree::node::ptr> >& list_tree) {
	tree::node::ptr new_p = valuation_functions::convert_tree_to_DNF(p);
	share_continuous_and_discrete_tree_recursion(new_p, list_tree);
}

}
