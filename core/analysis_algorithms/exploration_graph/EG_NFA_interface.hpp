/*
 * EG_NFA_interface.hpp
 *
 *  Created on: Apr 24, 2012
 *      Author: manish
 */

#ifndef EG_NFA_INTERFACE_HPP_
#define EG_NFA_INTERFACE_HPP_

#include "EG_NFA_interface.h"

namespace exploration_graph {

global_map_type EG_NFA_interface::init_global_map() {
	std::map<ExGraph*, local_abs_function_map_type> my_global_map;
	return my_global_map;
}

DFAutomaton EG_NFA_interface::create_DFA_from_path(const std::vector<
		ExGraphTransition> myTransitions) {

	DFAutomaton transVectorDFA;

	const AutomatonNode& newNode = transVectorDFA.createNode();
	const AutomatonNode* newNode_ptr = &newNode;

	transVectorDFA.setInitialNode(newNode);

	AutomatonNode* prevNode;
	prevNode = const_cast<AutomatonNode*> (newNode_ptr);

	for (std::vector<ExGraphTransition>::const_iterator tIt =
			myTransitions.begin(); tIt != myTransitions.end(); ++tIt) {

		std::string myNFALabel;
		ExGraphTransition_Info_ptr myInfo = (*tIt).getInfo();
		trans_type tType = myInfo->getType();

		if (tType == ExGraphTransition_Info::DISCRETE)
			myNFALabel = "1";
		else if (tType == ExGraphTransition_Info::TIME_ELAPSE)
			myNFALabel = "time";

		if (!myNFALabel.empty()) {
			const AutomatonNode& nextNode = transVectorDFA.createNode();
			const AutomatonNode* nextNode_ptr = &nextNode;
			transVectorDFA.addToAlphabet(myNFALabel);
			transVectorDFA.addEdge(*prevNode, nextNode, 1, myNFALabel);
			prevNode = const_cast<AutomatonNode*> (nextNode_ptr);
		}
	}
	transVectorDFA.setAcceptNode(*prevNode);
	return transVectorDFA;
}

/** Given a node, set of accepting nodes, an ExPath i.e., a list of transitions
 * , and the index in this path, this recursive function sets bad only those
 * transitions for this node, which lead to the bad path and also match the
 * transition, ExPath(index).
 *
 * O(m). */
bool check_and_set_bad_transitions(const ExNode_ptr exNode, const std::set<
		ExNode_ptr> acceptNodes, const std::vector<ExGraphTransition> ExPath,
		const std::vector<ExGraphTransition>::size_type index) {

	exNode->setStatus(ExGraphNode::IS_VISITED);

	bool is_bad = false;
	std::vector<ExGraphTransition> myTransitions = exNode->getTransitions();

	for (std::vector<ExGraphTransition>::iterator tIt = myTransitions.begin(); tIt
			!= myTransitions.end(); ++tIt) {

		bool local_is_bad = false;
		ExGraphTransition_Info_ptr myInfo = (*tIt).getInfo();
		ExNode_ptr succNode = access_exNode_method::access(
				(*tIt).getSuccessorNode());
		if (myInfo->getType() == ExGraphTransition_Info::CONTAINMENT) {
			if ((acceptNodes.count(succNode) > 0) && (index == ExPath.size()))
				local_is_bad = true;
			else if (succNode->getStatus() == ExGraphNode::IS_NOT_VISITED)
				local_is_bad = check_and_set_bad_transitions(succNode,
						acceptNodes, ExPath, index);
		} else {
			ExGraphTransition_Info_ptr info_at_index;
			ExGraphTransition ExPathTransition = ExPath.at(index);
			info_at_index = ExPathTransition.getInfo();
			if (info_at_index->IsEqual(myInfo)) {
				if ((acceptNodes.count(succNode) > 0) && ((index + 1)
						== ExPath.size())) {
					local_is_bad = true;
				} else if ((index + 1) < ExPath.size()) {
					local_is_bad = check_and_set_bad_transitions(succNode,
							acceptNodes, ExPath, (index + 1));
				} else
					local_is_bad = false;
			} else
				local_is_bad = false;
		}

		if (local_is_bad == true) {
			myInfo->setBadType(ExGraphTransition_Info::IS_BAD);
		}
		is_bad = is_bad || local_is_bad;
	}
	return is_bad;
}

NFAutomaton EG_NFA_interface::to_language(const ExGraph& exGraph,
		shortest_path_info_type EG_path_info) {

	std::vector < ExGraphTransition > ExPath = EG_path_info.Transitions;
	ExNode_ptr initNode_ExPath = EG_path_info.initialExNode;
	std::vector < ExNode_ptr > myExAcceptNodes = exGraph.getAcceptNodes();
	std::set < ExNode_ptr > myExAcceptNodesSet;

	for (std::vector<ExNode_ptr>::iterator aIt = myExAcceptNodes.begin(); aIt
			!= myExAcceptNodes.end(); ++aIt)
		myExAcceptNodesSet.insert(*aIt);

	std::vector<ExGraphTransition>::size_type myIndex = 0;
	bool is_bad = check_and_set_bad_transitions(initNode_ExPath,
			myExAcceptNodesSet, ExPath, myIndex);

	if(is_bad) {
		//!< do nothing
	}

	/** local abstraction function map for this exGraph. */
	local_abs_function_map_type local_abs_function_map;

	/** Create a new NFA. */
	NFAutomaton myNFA;

	/** ExNode to AutomatonNode mapping */
	std::map < ExNode_ptr, Node_ptr > ExNode_to_ANode_map;
	std::map<ExNode_ptr, Node_ptr>::iterator ExNode_to_ANode_mapIt;

	/** Retrieve the ExGraph information. */
	std::vector < ExNode_ptr > myExNodeList = exGraph.getNodeList();
	std::vector < ExNode_ptr > myExInitialNodes = exGraph.getInitialNodes();

	/** Create NFA nodes. */
	for (std::vector<ExNode_ptr>::iterator nIt = myExNodeList.begin(); nIt
			!= myExNodeList.end(); ++nIt) {
		(*nIt)->setStatus(ExGraphNode::IS_NOT_VISITED);
		const AutomatonNode& newNode_ref = myNFA.createNode();
		Node_ptr newNode = access_method::access(newNode_ref);
		ExNode_to_ANode_map.insert(std::pair<ExNode_ptr, Node_ptr>(*nIt,
				newNode));
	}

	/** Add Accepting nodes to NFA. */
	for (std::vector<ExNode_ptr>::iterator nIt = myExAcceptNodes.begin(); nIt
			!= myExAcceptNodes.end(); ++nIt) {
		ExNode_to_ANode_mapIt = ExNode_to_ANode_map.find(*nIt);
		if (ExNode_to_ANode_mapIt != ExNode_to_ANode_map.end()) {
			myNFA.setAcceptNode(*(ExNode_to_ANode_mapIt->second));
		}
	}

	std::set < ExNode_ptr > myExInitialNodeSet;

	/** Add initial nodes to NFA. */
	for (std::vector<ExNode_ptr>::iterator nIt = myExInitialNodes.begin(); nIt
			!= myExInitialNodes.end(); ++nIt) {
		ExNode_to_ANode_mapIt = ExNode_to_ANode_map.find(*nIt);
		if (ExNode_to_ANode_mapIt != ExNode_to_ANode_map.end()) {
			myNFA.setInitialNode(*(ExNode_to_ANode_mapIt->second));
			myExInitialNodeSet.insert(*nIt);
		}
	}

	/** for each initial node */
	while (!myExInitialNodeSet.empty()) {

		std::queue < ExNode_ptr > reachable_nodeQueue;

		ExNode_ptr initialNode;

		/** ExPath index for the ExPath initial node is set as 0. */
		if (myExInitialNodeSet.count(initNode_ExPath) > 0) {
			initialNode = initNode_ExPath;
		} else {
			initialNode = *(myExInitialNodeSet.begin());
		}

		if (initialNode->getStatus() == ExGraphNode::IS_NOT_VISITED) {
			initialNode->setStatus(ExGraphNode::IS_VISITED);
			reachable_nodeQueue.push(initialNode);
		}

		/** while the reachable node queue is not empty */
		while (!reachable_nodeQueue.empty()) {

			ExNode_ptr exNode = reachable_nodeQueue.front();

			reachable_nodeQueue.pop();

			ExNode_to_ANode_mapIt = ExNode_to_ANode_map.find(exNode);

			/** pointer to the starting node of an NFA edge. */
			Node_ptr sNode;
			if (ExNode_to_ANode_mapIt != ExNode_to_ANode_map.end())
				sNode = ExNode_to_ANode_mapIt->second;

			/** Get the transitions */
			std::vector < ExGraphTransition > myTransitions
					= exNode->getTransitions();

			ExTransition_to_label_container_type ExTransition_to_label_map;
			ExTransition_to_label_container_type::iterator
					ExTransition_to_label_mapIt;

			/** map to make sure that discrete transitions, for each node,
			 * with same ids map to same labels in an NFA. */
			std::map < transId_type, label_id_type
					> transitionId_to_labelId_map;
			std::map<transId_type, label_id_type>::iterator
					transitionId_to_labelId_mapIt;

			unsigned int myLabel = 2;
			/** for all transitions */
			for (std::vector<ExGraphTransition>::iterator tIt =
					myTransitions.begin(); tIt != myTransitions.end(); ++tIt) {

				/** Label for the corresponding NFA edge. */
				std::string myNFALabel;
				ExGraphTransition_Info_ptr myInfo = (*tIt).getInfo();
				ExNode_ptr mySuccNode = access_exNode_method::access(
						(*tIt).getSuccessorNode());
				trans_data_type tData = myInfo->getData();
				trans_type tType = myInfo->getType();
				trans_bad_type bad = myInfo->getBadType();

				label_id_type lId = label_id_type(-1);
				transId_type *tId = transId_type(0);

				if (bad == ExGraphTransition_Info::IS_BAD) {
					if (tType == ExGraphTransition_Info::DISCRETE) {
						tId = boost::get<transId_type>(&tData);
						myNFALabel = "1";
					} else if (tType == ExGraphTransition_Info::TIME_ELAPSE) {
						myNFALabel = "time_bad";
					} else
						myNFALabel = "EPS";
				}

				else if (bad == ExGraphTransition_Info::IS_NOT_BAD) {

					if (tType == ExGraphTransition_Info::DISCRETE) {
						tId = boost::get<transId_type>(&tData);
						transitionId_to_labelId_mapIt
								= transitionId_to_labelId_map.find(*tId);

						if (transitionId_to_labelId_mapIt
								== transitionId_to_labelId_map.end()) {
							std::stringstream ss;
							ss << myLabel;
							myNFALabel = ss.str();
							myLabel++;
						} else
							lId = transitionId_to_labelId_mapIt->second;
					} else if (tType == ExGraphTransition_Info::TIME_ELAPSE)
						myNFALabel = "time_not_bad";
					else
						myNFALabel = "EPS";
				}

				/** myNFALabel can be empty only when similar transition is already there for particular exNode. */
				if (!myNFALabel.empty())
					lId = label_cache::getLabelId(myNFALabel);

				/** Add label only if not already present in the label cache. */
				if (lId == -1)
					lId = label_cache::addLabel(myNFALabel);

				/** Only if a new NFALabel is computed for a discrete transition */
				if ((tType == ExGraphTransition_Info::DISCRETE)
						&& !myNFALabel.empty()) {
					transitionId_to_labelId_map.insert(std::pair<transId_type,
							label_id_type>(*tId, lId));
				}

				/** Add mapping. */
				ExTransition_to_label_map.insert(std::pair<
						ExGraphTransition_Info_ptr, label_id_type>(myInfo, lId));

				ExNode_to_ANode_mapIt = ExNode_to_ANode_map.find(mySuccNode);

				if (ExNode_to_ANode_mapIt != ExNode_to_ANode_map.end()) {
					Node_ptr succNode = ExNode_to_ANode_mapIt->second; // successor node of an edge

					myNFALabel = label_cache::getLabel(lId); // setting NFAlabel in case its not set above.

					myNFA.addToAlphabet(myNFALabel); // add label to NFA alphabet
					myNFA.addEdge(*sNode, *succNode, 1, myNFALabel); // add an edge

					/** If the node is already visited. Or say, all its transitions have already been
					 * checked, don't add this node to the queue. */
					if (mySuccNode->getStatus()
							== ExGraphNode::IS_NOT_VISITED) {
						reachable_nodeQueue.push(mySuccNode);
						mySuccNode->setStatus(ExGraphNode::IS_VISITED);
					}
				}
			}

			local_abs_function_map.insert(std::pair<ExNode_ptr,
					ExTransition_to_label_container_type>(exNode,
					ExTransition_to_label_map));

		} //!< while(!reachable_nodeQueue.empty())
		myExInitialNodeSet.erase(initialNode);
	} //!< for each initial node

	/** Reset status */
	for (std::vector<ExNode_ptr>::iterator nIt = myExNodeList.begin(); nIt
			!= myExNodeList.end(); ++nIt) {
		(*nIt)->setStatus(ExGraphNode::IS_NOT_VISITED);
	}

	ExGraph* exGraph_normal_ptr = const_cast<ExGraph*> (&exGraph);

	/** Add this exGraph map to the global map. */
	global_map.insert(std::pair<ExGraph*, local_abs_function_map_type>(
			exGraph_normal_ptr, local_abs_function_map));

	if (myNFA.is_empty())
		throw std::runtime_error("New encoded NFA is empty.");
	else
		return myNFA;
}

ExGraph EG_NFA_interface::from_language(const ExGraph& exGraph,
		const DFAutomaton& DFA) {
	NFAutomaton NFAOfDFA = DFA.convertToNFAStar();
	ExGraph EG;
	EG = from_language(exGraph, NFAOfDFA);
	return EG;
}

ExGraph EG_NFA_interface::from_language(const ExGraph& exGraph,
		const NFAutomaton& NFA) {

	/** label ids in an alphabet are mapped to their corresponding
	 * indices in the alphabet. */
	std::set < label_id_type > labelIds = NFA.getAlphabetIds();
	std::map<label_id_type, unsigned int> label_id_to_index_map;
	std::map<label_id_type, unsigned int>::iterator label_id_to_index_map_It;
	std::set<label_id_type>::size_type i = 0;
	for (std::set<label_id_type>::iterator lIt = labelIds.begin(); lIt
			!= labelIds.end(); ++lIt, ++i) {
		label_id_to_index_map.insert(std::pair<label_id_type, unsigned int>(
				*lIt, i));
	}

	/** NFAmap - maps each NFA node to a vector of nodes.
	 *
	 *  For a node, entry at Index i in its corresponding vector
	 *  represents the adjacent node reached from transition on
	 *  label at index i in NFAalphabet.
	 *
	 *  @note Useful because of the complete NFA.
	 */
	std::map < Node_ptr, std::vector<Node_ptr> > NFAmap;
	std::map<Node_ptr, std::vector<Node_ptr> >::iterator NFAmapIterator;

	std::vector < Node_ptr > NFANodeList = NFA.getNodeList();

	for (std::vector<Node_ptr>::iterator it = NFANodeList.begin(); it
			!= NFANodeList.end(); ++it) {

		/** Initialize NFAmap for each node with an emptyNFANodeVector, which
		 * will be filled later. */
		std::vector < Node_ptr > emptyNFANodeVector(labelIds.size());
		NFAmap.insert(std::pair<Node_ptr, std::vector<Node_ptr> >(*it,
				emptyNFANodeVector));

		NFAmapIterator = NFAmap.find(*it);

		std::vector < AutomatonEdge > myEdges = (*it)->getEdges();

		for (std::vector<AutomatonEdge>::iterator edgeIt = myEdges.begin(); edgeIt
				!= myEdges.end(); ++edgeIt) {

			Node_ptr succNode = access_method::access(
					edgeIt->getSuccessorNode());
			label_id_to_index_map_It = label_id_to_index_map.find(
					edgeIt->getEdgeLabelId());
			if (label_id_to_index_map_It != label_id_to_index_map.end())
				NFAmapIterator->second.at(label_id_to_index_map_It->second)
						= succNode;
		}
	}

	/** NFA has one initial state as it is formed from a DFA. */
	Node_ptr NFAinitialNode = NFA.getInitialNodes().at(0);
	std::set < Node_ptr > NFAacceptNodesSet;
	std::vector < Node_ptr > NFAacceptNodeList = NFA.getAcceptNodes();
	for (std::vector<Node_ptr>::iterator it = NFAacceptNodeList.begin(); it
			!= NFAacceptNodeList.end(); ++it)
		NFAacceptNodesSet.insert(*it);

	/** get initial node and accept nodes of the exploration graph. */
	std::vector < ExNode_ptr > ExinitialNodes = exGraph.getInitialNodes();
	std::set < ExNode_ptr > ExacceptNodesSet;
	std::vector < ExNode_ptr > ExacceptNodeList = exGraph.getAcceptNodes();
	for (std::vector<ExNode_ptr>::iterator it = ExacceptNodeList.begin(); it
			!= ExacceptNodeList.end(); ++it)
		ExacceptNodesSet.insert(*it);

	/** Iterator to the global map. */
	global_map_type::iterator global_mapIt;

	/** Get local_abs_function_map for this exploration graph. */
	local_abs_function_map_type local_abs_function_map;
	ExGraph* exGraph_normal_ptr = const_cast<ExGraph*> (&exGraph);
	global_mapIt = global_map.find(exGraph_normal_ptr);
	if (global_mapIt != global_map.end()) {
		local_abs_function_map = (*global_mapIt).second;
	}

	/** Create a new composite ExGraph */
	ExGraph compositeExGraph;

	/** Nodes reachable during parallel composition */
	std::queue < std::pair<ExNode_ptr, Node_ptr> > reachable_nodepairQueue;
	std::map < std::pair<ExNode_ptr, Node_ptr>, ExNode_ptr
			> NodeSet_to_Node_map;
	std::map<std::pair<ExNode_ptr, Node_ptr>, ExNode_ptr>::iterator
			NodeSet_to_Node_map_It1, NodeSet_to_Node_map_It2;

	std::vector<ExNode_ptr>::iterator ExinitialNodesIt = ExinitialNodes.begin();

	/** flags to check if either one or both of the nodes in a pair are accept nodes */
	bool acceptNode1, acceptNode2;
	acceptNode1 = false;
	acceptNode2 = false;

	/** Until all initial nodes have been traversed. */
	while (ExinitialNodesIt != ExinitialNodes.end()) {

		std::pair < ExNode_ptr, Node_ptr > reachable_nodepair;
		reachable_nodepair = std::pair<ExNode_ptr, Node_ptr>(*ExinitialNodesIt,
				NFAinitialNode);

		/** Get the iterator */
		NodeSet_to_Node_map_It1 = NodeSet_to_Node_map.find(reachable_nodepair);
		if (NodeSet_to_Node_map_It1 == NodeSet_to_Node_map.end()) {
			const ExGraphNode& newNode_ref = compositeExGraph.createNode();
			ExNode_ptr newNode = access_exNode_method::access(newNode_ref);
			compositeExGraph.setInitialNode(newNode_ref); // add initial node to exGraph.
			NodeSet_to_Node_map.insert(std::pair<
					std::pair<ExNode_ptr, Node_ptr>, ExNode_ptr>(
					reachable_nodepair, newNode));
			reachable_nodepairQueue.push(reachable_nodepair);
		} else
			compositeExGraph.setInitialNode(*(NodeSet_to_Node_map_It1->second));

		while (!reachable_nodepairQueue.empty()) {
			reachable_nodepair = reachable_nodepairQueue.front();
			ExNode_ptr firstNode;
			Node_ptr secondNode;
			firstNode = reachable_nodepair.first;
			secondNode = reachable_nodepair.second;
			reachable_nodepairQueue.pop();

			//!< node1 is an accepting node
			if (ExacceptNodesSet.count(firstNode) > 0)
				acceptNode1 = true;

			//!< node2 is an accepting node
			if (NFAacceptNodesSet.count(secondNode) > 0)
				acceptNode2 = true;

			/** Get the iterator */
			NodeSet_to_Node_map_It1 = NodeSet_to_Node_map.find(
					reachable_nodepair);

			if (NodeSet_to_Node_map_It1 == NodeSet_to_Node_map.end())
				throw std::runtime_error(
						"@Intersect: NodeSet is not there in the map.");

			if (acceptNode1 && acceptNode2)
				compositeExGraph.setAcceptNode(
						*(NodeSet_to_Node_map_It1->second));

			acceptNode1 = acceptNode2 = false;

			/** Get the vector for second node in the nodepair */
			NFAmapIterator = NFAmap.find(secondNode);

			local_abs_function_map_type::iterator local_abs_function_mapIt;
			local_abs_function_mapIt = local_abs_function_map.find(firstNode);
			ExTransition_to_label_container_type
					ExTransition_to_label_container;
			if (local_abs_function_mapIt != local_abs_function_map.end())
				ExTransition_to_label_container
						= local_abs_function_mapIt->second;

			std::vector < ExGraphTransition > myTransitions
					= firstNode->getTransitions();

			/** for each ExGraph transition */
			for (std::vector<ExGraphTransition>::iterator tIt =
					myTransitions.begin(); tIt != myTransitions.end(); ++tIt) {
				ExNode_ptr node1 = access_exNode_method::access(
						(*tIt).getSuccessorNode());
				ExGraphTransition_Info_ptr myInfo = (*tIt).getInfo();

				/** Get the labelId corresponding to this transition. */
				ExTransition_to_label_container_type::iterator
						ExTransition_to_label_containerIt;
				ExTransition_to_label_containerIt
						= ExTransition_to_label_container.find(myInfo);

				label_id_type myLabelId = label_id_type(-1);
				if (ExTransition_to_label_containerIt
						!= ExTransition_to_label_container.end())
					myLabelId = ExTransition_to_label_containerIt->second;

				std::set<label_id_type>::size_type l_index = 0;
				label_id_to_index_map_It
						= label_id_to_index_map.find(myLabelId);
				if (label_id_to_index_map_It != label_id_to_index_map.end()) {
					l_index = label_id_to_index_map_It->second;

					Node_ptr node2 = NFAmapIterator->second.at(l_index);

					if (node2 == NULL) {
						//!< std::cout << std::endl << "The node is null";
					} else {
						/** Create a new pair */
						reachable_nodepair = std::pair<ExNode_ptr, Node_ptr>(
								node1, node2);
						NodeSet_to_Node_map_It2 = NodeSet_to_Node_map.find(
								reachable_nodepair);

						/** If it is there in the map, add as a successor node */
						if (NodeSet_to_Node_map_It2
								!= NodeSet_to_Node_map.end()) {
							NodeSet_to_Node_map_It1->second->addTransition(
									*(NodeSet_to_Node_map_It2->second),
									myInfo->getType(), myInfo->getData());
						}//!< end-if

						else {

							/** else, create a new node corresponding to this pair
							 * add to the map and finally add as a successor node */
							const ExGraphNode& newNode_ref =
									compositeExGraph.createNode();
							ExNode_ptr newNode = access_exNode_method::access(
									newNode_ref);
							NodeSet_to_Node_map.insert(std::pair<std::pair<
									ExNode_ptr, Node_ptr>, ExNode_ptr>(
									reachable_nodepair, newNode));
							NodeSet_to_Node_map_It1->second->addTransition(
									*newNode, myInfo->getType(),
									myInfo->getData());
							reachable_nodepairQueue.push(reachable_nodepair);
						}//!< end-else if(NodeSet_to_Node_map_It2 != NodeSet_to_Node_map.end())
					}//!< end-else if (node2 == NULL)
				}//!< end if (label_id_to_index_map_It != label_id_to_index_map.end())
			}//!< end for each ExGraph transition
		}//!< end while (!reachable_nodepairQueue.empty())
		ExinitialNodesIt++;
	}//!< end while (ExinitialNodesIt != ExinitialNodes.end())

	ExGraph reverseEG = compositeExGraph.computeReverse(false);

	std::vector < ExNode_ptr > revExInitialNodes = reverseEG.getInitialNodes();
	std::queue < ExNode_ptr > revExReachNodes;
	for (std::vector<ExNode_ptr>::iterator initIt = revExInitialNodes.begin(); initIt
			!= revExInitialNodes.end(); ++initIt) {
		(*initIt)->setStatus(ExGraphNode::IS_VISITED);
		revExReachNodes.push(*initIt);
	}

	while (!revExReachNodes.empty()) {
		ExNode_ptr revExNode = revExReachNodes.front();
		revExReachNodes.pop();

		std::vector < ExGraphTransition > myTransitions
				= revExNode->getTransitions();
		for (std::vector<ExGraphTransition>::iterator tIt =
				myTransitions.begin(); tIt != myTransitions.end(); ++tIt) {
			const ExGraphNode& succNode = (*tIt).getSuccessorNode();
			ExNode_ptr succNode_shared_ptr = access_exNode_method::access(
					succNode);
			if (succNode_shared_ptr->getStatus()
					== ExGraphNode::IS_NOT_VISITED) {
				succNode_shared_ptr->setStatus(ExGraphNode::IS_VISITED);
				revExReachNodes.push(succNode_shared_ptr);
			}
		}
	}

	ExGraph finalEG = reverseEG.computeReverse(true);
	return finalEG;
}

}
#endif /* EG_NFA_INTERFACE_HPP_ */
