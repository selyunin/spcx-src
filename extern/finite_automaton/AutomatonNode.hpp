#ifndef AUTOMATONNODE_HPP_
#define AUTOMATONNODE_HPP_

#include "AutomatonNode.h"

namespace finite_automaton {

AutomatonNode::AutomatonNode() {
	myStatus = AutomatonNode::IS_NOT_VISITED;
}

AutomatonNode::~AutomatonNode() {
}

std::vector<AutomatonEdge> AutomatonNode::getEdges() const {
	return myEdges;
}

AutomatonNode::status AutomatonNode::getStatus() const {
	return myStatus;
}

void AutomatonNode::setStatus(const status st) {
	myStatus = st;
}
}

#endif /* AUTOMATONNODE_HPP_ */
