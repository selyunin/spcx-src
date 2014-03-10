/*
 * tree_traverser.cpp
 *
 *  Created on: Sep 7, 2009
 *      Author: frehse
 */

#include "utility/tree_traverser.h"

#include "utility/tree_node.h"

namespace tree {

tree_traverser::~tree_traverser() {}

void tree_traverser::traverse_preorder(node::ptr p)
	{
		if(eval(p))
		{
			if (binary_node* bp	= dynamic_cast<binary_node*>(p.get())) {
				traverse_preorder(bp->child1);
				traverse_preorder(bp->child2);
			}
		}
		else
			return;

	}

void tree_traverser::traverse_inorder(node::ptr p)
	{
		if (binary_node* bp
				= dynamic_cast<binary_node*>(p.get())) {
			traverse_inorder(bp->child1);
			eval(p);
			traverse_inorder(bp->child2);
		} else // it's a leaf node
		{
			eval(p);
		}
	}

void tree_traverser::traverse_postorder(node::ptr p)
	{
		if (binary_node* bp
				= dynamic_cast<binary_node*>(p.get())) {
			traverse_postorder(bp->child1);
			traverse_postorder(bp->child2);
			eval(p);
		} else // it's a leaf node
		{
			eval(p);
		}
	}


const_tree_traverser::~const_tree_traverser() {}

void const_tree_traverser::traverse_preorder(node::const_ptr p)
	{
		if(eval(p))
		{
			if (const binary_node* bp	= dynamic_cast<const binary_node*>(p.get())) {
				traverse_preorder(bp->child1);
				traverse_preorder(bp->child2);
			}
		}
		else
			return;

	}

void const_tree_traverser::traverse_inorder(node::const_ptr p)
	{
		if (const binary_node* bp
				= dynamic_cast<const binary_node*>(p.get())) {
			traverse_inorder(bp->child1);
			eval(p);
			traverse_inorder(bp->child2);
		} else // it's a leaf node
		{
			eval(p);
		}
	}

void const_tree_traverser::traverse_postorder(node::const_ptr p)
	{
		if (const binary_node* bp
				= dynamic_cast<const binary_node*>(p.get())) {
			traverse_postorder(bp->child1);
			traverse_postorder(bp->child2);
			eval(p);
		} else // it's a leaf node
		{
			eval(p);
		}
	}

}

