/*
 * convert_tree_to_DNF.cpp
 *
 *  Created on: Aug 27, 2009
 *      Author: gvincent
 */

#include "core/predicates/convert_tree_to_DNF.h"
#include "core/hybrid_automata/location_eq_node.h"

namespace valuation_functions {


bool is_DNF(const tree::node::ptr& p) {
	DNF_test g;
	g.traverse_preorder(p);
	return g.dnf;
}

bool is_OR_node(const tree::node::ptr& p) {
	if (boolean_node* b = dynamic_cast<boolean_node*>(p.get()))
					if (b->my_op == OR)
							return true;
	return false;
}

tree::node::ptr create_OR_AND_AND(const boolean_node* child1, const tree::node::ptr& child2) {
	tree::node::ptr child1_1 = convert_tree_to_DNF(child1->child1);
	tree::node::ptr child1_2 = convert_tree_to_DNF(child1->child2);

	tree::node::ptr new_child1 = tree::node::ptr(new boolean_node(AND,child1_1,child2));
	tree::node::ptr new_child2 = tree::node::ptr(new boolean_node(AND,child1_2,child2));
	tree::node::ptr new_node = tree::node::ptr(new boolean_node(OR,new_child1,new_child2));
	return new_node;
}

tree::node::ptr convert_tree_to_DNF(const tree::node::ptr& p) {
	return up_or_operator(down_not_operator(p));
}

tree::node::ptr up_or_operator(const tree::node::ptr& p) {

	tree::node::ptr new_node = p;
	if (boolean_node* q = dynamic_cast<boolean_node*>(p.get())) {
		tree::node::ptr child1 = up_or_operator(q->child1);
		tree::node::ptr child2 = up_or_operator(q->child2);
		new_node = tree::node::ptr(new boolean_node(q->my_op,child1,child2));

		if (q->my_op == AND) {
			if (is_OR_node(child1) && is_OR_node(child2)){
				boolean_node* c1 = dynamic_cast<boolean_node*>(child1.get());
				boolean_node* c2 = dynamic_cast<boolean_node*>(child2.get());
				new_node = tree::node::ptr(new boolean_node(OR,create_OR_AND_AND(c1, c2->child1),create_OR_AND_AND(c1, c2->child2)));
				new_node = up_or_operator(new_node);
			}
			else if (is_OR_node(child1)){
				boolean_node* c1 = dynamic_cast<boolean_node*>(child1.get());
				new_node = create_OR_AND_AND(c1, child2);
				new_node = up_or_operator(new_node);
			}
			else if (is_OR_node(child2)){
				boolean_node* c2 = dynamic_cast<boolean_node*>(child2.get());
				new_node = create_OR_AND_AND(c2, child1);
				new_node = up_or_operator(new_node);
			}
			/*
			 * no convert necessary if child1 and child2 are AND node
			 */
		}
	}

	return new_node;
}

tree::node::ptr down_not_operator(const tree::node::ptr& p) {

	tree::node::ptr new_node = p;
	if (boolean_node* q = dynamic_cast<boolean_node*>(p.get())) {
		if (q->my_op == NOT) {
			if (boolean_node* c = dynamic_cast<boolean_node*>(q->child1.get())) {
				tree::node::ptr child1 = down_not_operator(c->child1);
				if (c->my_op == NOT) {
					new_node = child1;
				}
				else {
					tree::node::ptr child2 = down_not_operator(c->child2);
					tree::node::ptr new_child1 = tree::node::ptr(new boolean_node(NOT,child1));
					tree::node::ptr new_child2 = tree::node::ptr(new boolean_node(NOT,child2));
					if (c->my_op == AND) {
						new_node = tree::node::ptr(new boolean_node(OR,new_child1,new_child2));
					}
					else if (c->my_op == OR) {
						new_node = tree::node::ptr(new boolean_node(AND,new_child1,new_child2));
					}
					else
						throw std::runtime_error("Unknown boolean operator");
				}
				new_node = down_not_operator(new_node);
			}
			else if (comparison_node* c = dynamic_cast<comparison_node*>(q->child1.get())){
				if(c->my_op == LT)
					new_node = tree::node::ptr(new comparison_node(GE,c->child1, c->child2));
				else if(c->my_op == GT)
					new_node = tree::node::ptr(new comparison_node(LE,c->child1, c->child2));
				else if(c->my_op == GE)
					new_node = tree::node::ptr(new comparison_node(LT,c->child1, c->child2));
				else if(c->my_op == LE)
					new_node = tree::node::ptr(new comparison_node(GT,c->child1, c->child2));
				else if(c->my_op == EQ){
					tree::node::ptr new_child1 = tree::node::ptr(new comparison_node(LT,c->child1, c->child2));
					tree::node::ptr new_child2 = tree::node::ptr(new comparison_node(GT,c->child1, c->child2));
					new_node = tree::node::ptr(new boolean_node(OR,new_child1,new_child2));
				}
				else
					throw std::runtime_error("Unknown comparison operator");
			}
			else if (hybrid_automata::location_eq_node* c = dynamic_cast<hybrid_automata::location_eq_node*>(q->child1.get())){
				new_node = tree::node::ptr(new hybrid_automata::location_eq_node(c->get_automaton_id(),c->get_location_id(),!(c->get_equal())));
			}
		}
		else {
			tree::node::ptr child1 = down_not_operator(q->child1);
			tree::node::ptr child2 = down_not_operator(q->child2);
			new_node = tree::node::ptr(new boolean_node(q->my_op,child1,child2));
		}
	}
	return new_node;
}

}
