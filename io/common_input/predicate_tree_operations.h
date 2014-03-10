/*
 * predicate_tree_operations.h
 *
 *  Created on: Sep 4, 2009
 *      Author: frehse
 */

#ifndef PREDICATE_TREE_OPERATIONS_H_
#define PREDICATE_TREE_OPERATIONS_H_

#include "boost/shared_ptr.hpp"

/** Forward declarations */
namespace tree {
class node;
typedef boost::shared_ptr<node> node_ptr;
}

std::pair<tree::node_ptr, tree::node_ptr> divide_tree(tree::node_ptr p);

#endif /* PREDICATE_TREE_OPERATIONS_H_ */
