/*
 * tree_node_cache.h
 *
 *  Created on: Jun 28, 2011
 *      Author: goyal
 */

#ifndef TREE_NODE_CACHE_H_
#define TREE_NODE_CACHE_H_

#include <string>
#include <map>
#include <vector>
#include <stdexcept>
#include <iostream>
#include "boost/shared_ptr.hpp"
#include "../utility/tree_node_id.h"

using namespace std;
/** Forward declarations */
namespace tree {
class node;
typedef boost::shared_ptr<node> node_ptr;
//typedef boost::shared_ptr<const discrete_set> const_ptr;
}

namespace tree {

/** Stores automata globally.
 */

template<typename scalar_type>
class tree_node_cache {

public:
	/** Add a node to the cache. */

	static void add_node(node_ptr n, const scalar_type val) {
		std::cout << "\nI'm in add node";
		tree_node_id id;
		id = get_new_id();
		tree_node_id_ptr_map[id] = n;
		tree_node_id_val_map[id] = val;
	}

	static tree_node_id has_node(const node_ptr& n) {
		for (container_type::const_iterator it = tree_node_id_ptr_map.begin(); it
					!= tree_node_id_ptr_map.end(); ++it) {
				if (it->second == n)
					return it->first;
		}
		return 0;
	}

	static scalar_type get_val(const tree_node_id id) {
		std::cout << "\nI'm in get val";
		return tree_node_id_val_map[id];

	}

	/** Output as a stream of characters.
	 */
	static void print(std::ostream& os);


	/** Retrieve the set of all node ids in the cache. */
	static std::set<tree_node_id> get_nodes();

	static tree_node_id get_new_id() {
		return ++highest_id;
	}


private:
	static tree_node_id highest_id;
	typedef std::map<tree_node_id, node_ptr> container_type;
	static container_type tree_node_id_ptr_map;
	typedef std::map<tree_node_id, scalar_type> container_type1;
	static container_type1 tree_node_id_val_map;
};

}


#endif /* TREE_NODE_CACHE_H_ */
