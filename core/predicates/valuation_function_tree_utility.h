#ifndef VALUATION_FUNCTION_TREE_UTILITY_H_
#define VALUATION_FUNCTION_TREE_UTILITY_H_

#include "math/vdom/variable.h"
//#include "../utility/tree_node.h" // needed for tree_traverser (should be changed)
//#include "valuation_function_tree.h"
#include "valuation_function_tree_nodes.h"

#include <boost/shared_ptr.hpp>

/** Forward declarations */
namespace tree {
class node;
typedef boost::shared_ptr<node> node_ptr;
typedef boost::shared_ptr<const node> node_const_ptr;
}

namespace valuation_functions {

/** Returns the boolean and of two predicates
 *
 * If one of the predicates is null, return the other one. */
tree::node_ptr boolean_and(const tree::node_ptr& p1, const tree::node_ptr& p2);

/** Return the ids of all variables in the predicate.
 *
 * The ids are given as is, without any modification (primedness etc.).*/
variable_id_set get_variable_ids(tree::node_const_ptr p);

/** Get variable_ids irrespective of their prime count and
 * return their unprimed ids.
 */
variable_id_set get_unprimed_variable_ids(tree::node_ptr p);

/** Get variable_ids of primedness prime_count and
 * return their unprimed ids. */
variable_id_set get_primed_variable_ids(tree::node_ptr p, unsigned int prime_count);

/** Increase the primedness of variables in the tree by delta.
 * The default value for delta is 1. */
void increase_primedness(tree::node_ptr p, unsigned int delta=1);

/** Add const constraints to a transform predicate
 *
 * The constraint x' == x is added for every variable x in vis.
 * Note that the variable in vis is not primed.
 */
void add_const_relation_constraints(tree::node_ptr& p, const variable_id_set& vis);

/** Add const constraints to a dynamics predicate
 *
 * The constraint x' == 0 is added for every variable x in vis.
 * Note that the variable in vis is not primed.
 */
void add_const_dynamics_constraints(tree::node_ptr& p, const variable_id_set& vis);

}

#endif /*VALUATION_FUNCTION_TREE_UTILITY_H_*/
