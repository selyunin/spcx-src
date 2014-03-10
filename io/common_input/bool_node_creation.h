/*
 * bool_node_creation.h
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#ifndef BOOL_NODE_CREATION_H_
#define BOOL_NODE_CREATION_H_

#include "boost/shared_ptr.hpp"

/** Forward declarations */
namespace tree {
class node;
typedef boost::shared_ptr<node> node_ptr;
}

namespace predicate_parser {

/** Create a const_node from string, taking the type from parse_type_chooser */
tree::node_ptr create_bool_const_node(std::string s);

/** test if p is a bool const_node*/
bool is_bool_const_node(const tree::node_ptr& p);

/** Return the value of a bool const_node
 *
 * Throws if unsuccessful.
 */
bool get_bool_const_node_value(const tree::node_ptr& p);

}

#endif /* BOOL_NODE_CREATION_H_ */
