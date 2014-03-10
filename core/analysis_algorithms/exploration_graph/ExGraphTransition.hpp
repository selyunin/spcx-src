/*
 * ExGraphTransition.hpp
 *
 *  Created on: Apr 17, 2012
 *      Author: manish
 */

#ifndef EXGRAPHTRANSITION_HPP_
#define EXGRAPHTRANSITION_HPP_

/** namespace contains a method for converting an ExNode
 * reference to its shared pointer.*/
namespace access_exNode_method {
using namespace exploration_graph;
typedef boost::shared_ptr<ExGraphNode> ExNode_ptr;
ExNode_ptr access(const ExGraphNode& node) {
	ExGraphNode* node_ptr = const_cast<ExGraphNode*>(&node);
	ExNode_ptr node_shared_ptr = node_ptr->get_ptr();
	return node_shared_ptr;
}
}
;

#include "ExGraphTransition.h"

namespace exploration_graph {

ExGraphTransition::ExGraphTransition() {
}

ExGraphTransition::~ExGraphTransition() {
}

ExGraphTransition::ExGraphTransition(const ExGraphNode& succNode,
		const trans_type& tType, const trans_data_type& tData) {

	ExGraphTransition_Info_ptr transInfo(new ExGraphTransition_Info(tType, tData));
	successorNode = access_exNode_method::access(succNode);
	myInfo = transInfo;
}
}
#endif /* EXGRAPHTRANSITION_HPP_ */
