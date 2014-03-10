/*
 * node_print_visitor.h
 *
 *  Created on: Aug 28, 2009
 *      Author: gvincent
 */

#ifndef NODE_PRINT_VISITOR_H_
#define NODE_PRINT_VISITOR_H_

#include "core/predicates/valuation_function_tree_nodes.h"

namespace valuation_functions {

std::string operator_to_string(const arithmetic_node* p);
std::string operator_to_string(const boolean_node* p);
std::string operator_to_string(const comparison_node* p);

void print_node(std::ostream& os, const tree::node* p);

}

std::ostream& operator<<(std::ostream& os, const tree::node& p);
std::ostream& operator<<(std::ostream& os, const tree::node::ptr& p);
std::ostream& print_as_predicate(std::ostream& os, const tree::node::ptr& p);
std::ostream& print_as_arithmetic(std::ostream& os, const tree::node::ptr& p);

#endif /* NODE_PRINT_VISITOR_H_ */
