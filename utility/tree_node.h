#ifndef TREE_NODE_H_
#define TREE_NODE_H_

#include <boost/shared_ptr.hpp>
//#include "../utility/tree_node_id.h"
//#include "../utility/shared_ptr_user.h"
//#include "../valuation_function/node_print_visitor.h"

namespace tree {
class node_visitor;

//class node : public shared_ptr_user<node> {
class node {
public:
	typedef boost::shared_ptr<node> ptr;
	typedef boost::shared_ptr<const node> const_ptr;

	node() {
	}
	;
	virtual ~node() {
	}
	;

	virtual void accept (tree::node_visitor& v) const = 0;

	//virtual const tree_node_id& get_id() const;

	static node::ptr null_node() { return boost::shared_ptr<node>(); };


};

class binary_node : public node {
public:
	binary_node(const node::ptr& child_one, const node::ptr& child_two) :
		child1(child_one), child2(child_two) {

	}
	;
	virtual ~binary_node() {
	}
	;
	node::ptr child1;
	node::ptr child2;

};

}

#include "utility/tree_node_declarations.h"
#include "utility/tree_node_visitor.h"



#endif /*TREE_NODE_H_*/
