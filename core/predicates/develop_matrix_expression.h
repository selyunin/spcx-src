/*
 * develop_matrix_expression.h
 *
 *  Created on: Oct 12, 2010
 *      Author: gvincent
 */

#ifndef DEVELOP_MATRIX_EXPRESSION_H_
#define DEVELOP_MATRIX_EXPRESSION_H_

#include "boost/shared_ptr.hpp"

/** Forward declarations */
namespace tree {
class node;
typedef boost::shared_ptr<node> node_ptr;
}

namespace valuation_functions {

/**
 * Develop matrix assignment in a tree.
 * E. g. "x' = A * x + b" (with x(2), A(2,2) and b(2) become:
 * "x(1)' = A(1, 1) * x(1) + A(1, 2) * x(2) + b(1) && x(2)' = A(2, 1) * x(1) + A(2, 2) * x(2) + b(2)".
 * Other nodes do not change.
 */
tree::node_ptr develop_assignment_matrix_tree(const tree::node_ptr& p);
}

#endif /* DEVELOP_MATRIX_EXPRESSION_H_ */
