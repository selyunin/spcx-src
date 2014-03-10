/*
 * share_discrete_and_continuous_tree.h
 *
 *  Created on: Aug 27, 2009
 *      Author: gvincent
 */

#ifndef SHARE_DISCRETE_AND_CONTINUOUS_TREE_H_
#define SHARE_DISCRETE_AND_CONTINUOUS_TREE_H_

#include <vector>

/** Forward declarations */
namespace tree {
class node;
typedef boost::shared_ptr<node> node_ptr;
}

namespace hybrid_automata {

void share_continuous_and_discrete_tree(const tree::node_ptr& p,
		std::vector<std::pair<tree::node_ptr, tree::node_ptr> >& list_tree);

}

#endif /* SHARE_DISCRETE_AND_CONTINUOUS_TREE_H_ */
