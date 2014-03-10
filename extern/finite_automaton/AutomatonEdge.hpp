#ifndef AUTOMATONEDGE_HPP_
#define AUTOMATONEDGE_HPP_

/** namespace contains a method for converting an automaton-node
 * reference to its shared pointer.*/
namespace access_method {
using namespace finite_automaton;
typedef boost::shared_ptr<AutomatonNode> Node_ptr;
Node_ptr access(const AutomatonNode& node) {
	AutomatonNode* node_ptr = const_cast<AutomatonNode*> (&node);
	Node_ptr node_shared_ptr = node_ptr->get_ptr();
	return node_shared_ptr;
}
}
;

#include "AutomatonEdge.h"

namespace finite_automaton {

AutomatonEdge::AutomatonEdge() {
}

AutomatonEdge::AutomatonEdge(const AutomatonNode& eSuccNode,
		const cost_type& eCost, const label_id_type& eLabelId) {

	AutomatonNode* node_ptr = const_cast<AutomatonNode*> (&eSuccNode);
	Node_ptr node_shared_ptr = node_ptr->get_ptr();

	successorNode = node_shared_ptr;
	edgeCost = eCost;
	edgeLabel_id = eLabelId;
}

AutomatonEdge::~AutomatonEdge() {
}

const AutomatonNode& AutomatonEdge::getSuccessorNode() const {
	return *successorNode;
}

const AutomatonEdge::cost_type& AutomatonEdge::getEdgeCost() const {
	return edgeCost;
}

const std::string AutomatonEdge::getEdgeLabel() const {
	return label_cache::getLabel(edgeLabel_id);
}
}
#endif /* AUTOMATONEDGE_HPP_ */
