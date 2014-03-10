/*
 * scalar_node_creation.h
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#ifndef SCALAR_NODE_CREATION_H_
#define SCALAR_NODE_CREATION_H_

#include "boost/shared_ptr.hpp"

/** Forward declarations */
namespace tree {
class node;
typedef boost::shared_ptr<node> node_ptr;
}

namespace parser {
class symbol;
}


namespace predicate_parser {

/** Create a const_node from string, taking the type from parse_type_chooser */
tree::node_ptr create_scalar_const_node(const std::string& s);

/** Create a const_node from a matrix of any, (only Rational type is implemented for matrix const_node) */
tree::node_ptr create_matrix_scalar_const_node(const parser::symbol& s);

/** test if p is a scalar const_node*/
bool is_scalar_const_node(const tree::node_ptr& p);
};

#endif /* SCALAR_NODE_CREATION_H_ */
