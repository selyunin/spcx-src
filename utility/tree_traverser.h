/*
 * tree_traverser.h
 *
 *  Created on: Sep 7, 2009
 *      Author: frehse
 */

#ifndef TREE_TRAVERSER_H_
#define TREE_TRAVERSER_H_

#include <boost/shared_ptr.hpp>

/** Forward declarations */
namespace tree {
class node;
typedef boost::shared_ptr<node> node_ptr;
typedef boost::shared_ptr<const node> node_const_ptr;
}

namespace tree {

class tree_traverser {
public:
	virtual ~tree_traverser();

	virtual bool eval(node_ptr p) = 0;

	/** Traverse the tree, children last. */
	virtual void traverse_preorder(node_ptr p);

	/** Traverse the tree, left child first, then *this, then right child. */
	virtual void traverse_inorder(node_ptr p);

	/** Traverse the tree, children first. */
	virtual void traverse_postorder(node_ptr p);
};

class const_tree_traverser {
public:
	virtual ~const_tree_traverser();

	virtual bool eval(node_const_ptr p) = 0;

	/** Traverse the tree, children last. */
	virtual void traverse_preorder(node_const_ptr p);

	/** Traverse the tree, left child first, then *this, then right child. */
	virtual void traverse_inorder(node_const_ptr p);

	/** Traverse the tree, children first. */
	virtual void traverse_postorder(node_const_ptr p);
};

}

#endif /* TREE_TRAVERSER_H_ */
