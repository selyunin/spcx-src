/*
 * ExGraphNode.hpp
 *
 *  Created on: Apr 17, 2012
 *      Author: manish
 */

#ifndef EXGRAPHNODE_HPP_
#define EXGRAPHNODE_HPP_

#include "ExGraphNode.h"

namespace exploration_graph {

ExGraphNode::ExGraphNode() {
	myStatus = ExGraphNode::IS_NOT_VISITED;
}

ExGraphNode::~ExGraphNode(void) {
}

std::vector<ExGraphTransition> ExGraphNode::getTransitions() const {
	return myTransitions;
}

ExGraphNode::status ExGraphNode::getStatus() const {
	return myStatus;
}

void ExGraphNode::setStatus(const status st) {
	myStatus = st;
}

void ExGraphNode::addTransition(const ExGraphNode& succNode, const trans_type& tType,
		const trans_data_type& tData) {

	ExGraphTransition trans(succNode, tType, tData);
	myTransitions.push_back(trans);
}
}
#endif /* EXGRAPHNODE_HPP_ */
