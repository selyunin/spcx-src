#ifndef NFAUTOMATON_HPP_
#define NFAUTOMATON_HPP_

#include "Theset.h"
#include "NFAutomaton.h"

namespace finite_automaton {

NFAutomaton::NFAutomaton() {
}

NFAutomaton::~NFAutomaton() {
}

NFAutomaton::NFAutomaton(const NFAutomaton& automaton) {

	std::vector < Node_ptr > NodeList = automaton.getNodeList();

	std::set < label_cache::label_id_type > alphabetIds
			= automaton.getAlphabetIds();

	this->setAlphabetIds(alphabetIds);
	std::map < Node_ptr, Node_ptr > Nodes_map;
	std::map<Node_ptr, Node_ptr>::iterator Nodes_mapIt;
	for (std::vector<Node_ptr>::const_iterator nIt = NodeList.begin(); nIt
			!= NodeList.end(); ++nIt) {
		const AutomatonNode& newNodeObj_ref = this->createNode();
		Node_ptr newNode = access_method::access(newNodeObj_ref);
		Nodes_map.insert(std::pair<Node_ptr, Node_ptr>(*nIt, newNode));
	}

	for (std::vector<Node_ptr>::const_iterator nIt = NodeList.begin(); nIt
			!= NodeList.end(); ++nIt) {
		Node_ptr sNode;
		Nodes_mapIt = Nodes_map.find(*nIt);
		if (Nodes_mapIt != Nodes_map.end())
			sNode = Nodes_mapIt->second;
		std::vector < AutomatonEdge > myEdges = (*nIt)->getEdges();
		for (std::vector<AutomatonEdge>::iterator eIt = myEdges.begin(); eIt
				!= myEdges.end(); ++eIt) {

			const AutomatonNode& succNodeObj_ref = (*eIt).getSuccessorNode();
			Node_ptr succNode = access_method::access(succNodeObj_ref);
			cost_type cost = (*eIt).getEdgeCost();
			label_cache::label_id_type lId = (*eIt).getEdgeLabelId();

			Node_ptr tNode;

			Nodes_mapIt = Nodes_map.find(succNode);
			if (Nodes_mapIt != Nodes_map.end())
				tNode = Nodes_mapIt->second;
			sNode->addSuccNode(*tNode, cost, lId);
		}
	}
	std::vector < Node_ptr > acceptNodes = automaton.getAcceptNodes();
	for (std::vector<Node_ptr>::const_iterator nIt = acceptNodes.begin(); nIt
			!= acceptNodes.end(); ++nIt) {
		Nodes_mapIt = Nodes_map.find(*nIt);
		if (Nodes_mapIt != Nodes_map.end())
			setAcceptNode(*(Nodes_mapIt->second));
	}
	std::vector < Node_ptr > initialNodes = automaton.getInitialNodes();
	for (std::vector<Node_ptr>::const_iterator nIt = initialNodes.begin(); nIt
			!= initialNodes.end(); ++nIt) {
		Nodes_mapIt = Nodes_map.find(*nIt);
		if (Nodes_mapIt != Nodes_map.end())
			this->setInitialNode(*(Nodes_mapIt->second));
	}
}

NFAutomaton NFAutomaton::concatenate(const NFAutomaton& secondNFA) const {

	NFAutomaton firstNFA_copy(*this);
	NFAutomaton secondNFA_copy(secondNFA);

	/** Create a new NFA */
	NFAutomaton concatenated_NFA;

	/** Set alphabet of concatenated NFA */
	std::set<label_cache::label_id_type> labelIdSet, secondNFAlabelIdSet;

	labelIdSet = firstNFA_copy.getAlphabetIds();

	/** "EPS" label has to be there in concatenated NFA alphabet
	 * because NULL transitions were added to it above. */
	labelIdSet.insert(label_cache::label_id_type(0));

	secondNFAlabelIdSet = secondNFA_copy.getAlphabetIds();

	labelIdSet.insert(secondNFAlabelIdSet.begin(), secondNFAlabelIdSet.end());

	concatenated_NFA.setAlphabetIds(labelIdSet);

	/** Merge multiple initial nodes, if any */
	firstNFA_copy.mergeInitialNodes(true);
	secondNFA_copy.mergeInitialNodes(true);

	std::vector < Node_ptr > myAcceptNodes = firstNFA_copy.getAcceptNodes();
	std::set < Node_ptr > myAcceptNodesSet;
	for (std::vector<Node_ptr>::iterator aIt = myAcceptNodes.begin(); aIt
			!= myAcceptNodes.end(); ++aIt)
		myAcceptNodesSet.insert(*aIt);

	std::vector < Node_ptr > myInitialNodes = firstNFA_copy.getInitialNodes();

	std::vector < Node_ptr > myNodes = firstNFA_copy.getNodeList();

	std::map < Node_ptr, Node_ptr > origNodetoNFANodemap;
	std::map<Node_ptr, Node_ptr>::iterator origNodetoNFANodemapIt;

	std::vector < Node_ptr > temp_acceptNodes;
	std::vector < Node_ptr > newNFANodes;

	/** Add nodes from first NFA to newly created concatenated_NFA.
	 *
	 * @note Initial node of first NFA corresponds to the initial node
	 * of new concatenated NFA. Whereas, firstNFA accepting nodes are
	 * stored in a temporary vector from which, NULL transitions
	 * will be added to second NFA's initial node.
	 */

	for (std::vector<Node_ptr>::size_type i = 0; i < myNodes.size(); ++i) {
		const AutomatonNode& newNodeObj_ref = concatenated_NFA.createNode();
		Node_ptr newNode = access_method::access(newNodeObj_ref);
		newNFANodes.push_back(newNode);
		if (myAcceptNodesSet.count(myNodes[i]) > 0)
			temp_acceptNodes.push_back(newNode);
		if (myInitialNodes[0] == myNodes[i])
			concatenated_NFA.setInitialNode(*newNode);
		origNodetoNFANodemap.insert(std::pair<Node_ptr, Node_ptr>(myNodes[i],
				newNode));
	}

	/** Copy transitions from first NFA */
	for (std::vector<Node_ptr>::size_type i = 0; i < myNodes.size(); ++i) {
		std::vector < AutomatonEdge > myEdges = myNodes[i]->getEdges();

		for (std::vector<AutomatonEdge>::size_type j = 0; j < myEdges.size(); ++j) {
			const AutomatonNode& succNodeObj_ref =
					myEdges[j].getSuccessorNode();
			Node_ptr succNode = access_method::access(succNodeObj_ref);
			origNodetoNFANodemapIt = origNodetoNFANodemap.find(succNode);

			if (origNodetoNFANodemapIt != origNodetoNFANodemap.end())
				newNFANodes[i]->addSuccNode(*(origNodetoNFANodemapIt->second),
						myEdges[j].getEdgeCost(), myEdges[j].getEdgeLabelId());
		}
	}

	/** Clear data structures to be reused for second NFA. */
	myAcceptNodes.clear();
	myAcceptNodesSet.clear();
	origNodetoNFANodemap.clear();
	myInitialNodes.clear();
	myNodes.clear();
	newNFANodes.clear();

	/** Second NFA accept Nodes */
	myAcceptNodes = secondNFA_copy.getAcceptNodes();

	for (std::vector<Node_ptr>::iterator aIt = myAcceptNodes.begin(); aIt
			!= myAcceptNodes.end(); ++aIt)
		myAcceptNodesSet.insert(*aIt);

	myInitialNodes = secondNFA_copy.getInitialNodes();

	/** Initial node is made temporary */
	Node_ptr temp_initNode;

	myNodes = secondNFA_copy.getNodeList();

	/** Add nodes from second NFA to the new NFAutomaton.
	 *
	 * @note accepting nodes from second NFA corresponds to the accepting nodes
	 * of new concatenated NFA. Whereas, secondNFA initial node is stored in a
	 * temp variable, temp_initNode.
	 */
	for (std::vector<Node_ptr>::size_type i = 0; i < myNodes.size(); ++i) {
		const AutomatonNode& newNodeObj_ref = concatenated_NFA.createNode();
		Node_ptr newNode = access_method::access(newNodeObj_ref);
		newNFANodes.push_back(newNode);
		if (myAcceptNodesSet.count(myNodes[i]) > 0)
			concatenated_NFA.setAcceptNode(*newNode);
		if (myInitialNodes[0] == myNodes[i])
			temp_initNode = newNode;
		origNodetoNFANodemap.insert(std::pair<Node_ptr, Node_ptr>(myNodes[i],
				newNode));
	}

	/** Copy transitions from secondNFA */
	for (std::vector<Node_ptr>::size_type i = 0; i < myNodes.size(); ++i) {
		std::vector < AutomatonEdge > myEdges = myNodes[i]->getEdges();

		for (std::vector<AutomatonEdge>::size_type j = 0; j < myEdges.size(); ++j) {
			const AutomatonNode& succNodeObj_ref =
					myEdges[j].getSuccessorNode();
			Node_ptr succNode = access_method::access(succNodeObj_ref);
			origNodetoNFANodemapIt = origNodetoNFANodemap.find(succNode);
			if (origNodetoNFANodemapIt != origNodetoNFANodemap.end())
				newNFANodes[i]->addSuccNode(*(origNodetoNFANodemapIt->second),
						myEdges[j].getEdgeCost(), myEdges[j].getEdgeLabelId());
		}
	}

	/** Add NULL transitions from firstNFA accepting nodes to secondNFA initial node. */
	for (std::vector<Node_ptr>::iterator it = temp_acceptNodes.begin(); it
			!= temp_acceptNodes.end(); ++it)
		(*it)->addSuccNode(*temp_initNode, 1, label_cache::label_id_type(0));

	return concatenated_NFA;

}

NFAutomaton NFAutomaton::computeKleeneStar() const {

	NFAutomaton thisNFA_copy(*this);

	thisNFA_copy.mergeInitialNodes(true);

	std::vector < Node_ptr > myInitialNodes;

	myInitialNodes = thisNFA_copy.getInitialNodes();

	NFAutomaton KleeneOfNFA;

	std::set < label_cache::label_id_type > mylabelIds
			= thisNFA_copy.getAlphabetIds();
	KleeneOfNFA.setAlphabetIds(mylabelIds);

	/** Only if EPS transition is not already there */
	if (!(mylabelIds.count(label_cache::label_id_type(0)) > 0))
		KleeneOfNFA.addToAlphabet("EPS");

	std::vector < Node_ptr > myAcceptNodes = thisNFA_copy.getAcceptNodes();
	std::set < Node_ptr > myAcceptNodesSet;
	for (std::vector<Node_ptr>::iterator aIt = myAcceptNodes.begin(); aIt
			!= myAcceptNodes.end(); ++aIt)
		myAcceptNodesSet.insert(*aIt);

	std::vector < Node_ptr > myNodes = thisNFA_copy.getNodeList();

	std::map < Node_ptr, Node_ptr > origNodetoKleeneNodemap;
	std::map<Node_ptr, Node_ptr>::iterator origNodetoKleeneNodemapIt;

	Node_ptr oldInitialNode;

	for (std::vector<Node_ptr>::size_type i = 0; i < myNodes.size(); ++i) {
		const AutomatonNode& newNodeObj_ref = KleeneOfNFA.createNode();
		Node_ptr newNode = access_method::access(newNodeObj_ref);
		if (myAcceptNodesSet.count(myNodes[i]) > 0)
			KleeneOfNFA.setAcceptNode(*newNode);
		if (myInitialNodes[0] == myNodes[i])
			oldInitialNode = newNode;
		origNodetoKleeneNodemap.insert(std::pair<Node_ptr, Node_ptr>(
				myNodes[i], newNode));
	}

	std::vector < Node_ptr > NFANodes = KleeneOfNFA.getNodeList();

	for (std::vector<Node_ptr>::size_type i = 0; i < myNodes.size(); ++i) {
		std::vector < AutomatonEdge > myEdges = myNodes[i]->getEdges();

		for (std::vector<AutomatonEdge>::size_type j = 0; j < myEdges.size(); ++j) {
			const AutomatonNode& succNodeObj_ref =
					myEdges[j].getSuccessorNode();
			Node_ptr succNode = access_method::access(succNodeObj_ref);
			origNodetoKleeneNodemapIt = origNodetoKleeneNodemap.find(succNode);

			if (origNodetoKleeneNodemapIt != origNodetoKleeneNodemap.end())
				NFANodes[i]->addSuccNode(*(origNodetoKleeneNodemapIt->second),
						myEdges[j].getEdgeCost(), myEdges[j].getEdgeLabelId());
		}
	}

	const AutomatonNode& newNodeObj_ref = KleeneOfNFA.createNode();
	Node_ptr newNode = access_method::access(newNodeObj_ref);

	/** It has only one initial state because multiple initial states (if any)
	 * have already been merged in the beginning. */
	newNode->addSuccNode(*oldInitialNode, 1, label_cache::label_id_type(0));

	myAcceptNodes.clear();
	myAcceptNodes = KleeneOfNFA.getAcceptNodes();

	for (std::vector<Node_ptr>::size_type i = 0; i < myAcceptNodes.size(); ++i)
		myAcceptNodes[i]->addSuccNode(*oldInitialNode, 1,
				label_cache::label_id_type(0));

	/** make newnode as initial node */
	KleeneOfNFA.setInitialNode(*newNode);

	/** Add newnode to the AcceptNodes' list */
	KleeneOfNFA.setAcceptNode(*newNode);

	return KleeneOfNFA;

}

DFAutomaton NFAutomaton::convertToDFA() const {

	NFAutomaton thisNFA_copy(*this);

	thisNFA_copy.mergeInitialNodes(true);

	std::vector < Node_ptr > NodeList = thisNFA_copy.getNodeList();

	/** For each node in an NFA, each transition label is mapped to a set */
	std::map < Node_ptr, std::map<label_cache::label_id_type,
			std::set<Node_ptr> > > globalmap;
	std::map<Node_ptr,
			std::map<label_cache::label_id_type, std::set<Node_ptr> > >::iterator
			globalmapIterator;
	std::map<label_cache::label_id_type, std::set<Node_ptr> >::iterator
			mymapIterator;

	for (std::vector<Node_ptr>::size_type i = 0; i < NodeList.size(); ++i) {

		std::vector < AutomatonEdge > myEdges = NodeList[i]->getEdges();

		/** Each transition label is mapped to a set */
		std::map < label_cache::label_id_type, std::set<Node_ptr> > mymap;

		for (std::vector<AutomatonEdge>::size_type j = 0; j < myEdges.size(); ++j) {

			label_cache::label_id_type lId = myEdges[j].getEdgeLabelId();

			mymapIterator = mymap.find(lId);

			const AutomatonNode& succNodeObj_ref =
					myEdges[j].getSuccessorNode();
			Node_ptr succNode = access_method::access(succNodeObj_ref);

			std::set < Node_ptr > aNodeSet;

			aNodeSet.insert(succNode);

			if (mymapIterator == mymap.end() && !aNodeSet.empty()) {
				mymap.insert(std::pair<label_cache::label_id_type, std::set<
						Node_ptr> >(lId, aNodeSet));
			}

			else
				(*mymapIterator).second.insert(aNodeSet.begin(), aNodeSet.end());

		}
		globalmap.insert(std::pair<Node_ptr, std::map<
				label_cache::label_id_type, std::set<Node_ptr> > >(NodeList[i],
				mymap));
	}

	/** Compute Null closures
	 *
	 * Since, in the beginning we have merged multiple initial states, we have just
	 * one NFA initial state now. */
	Node_ptr iNode = thisNFA_copy.getInitialNodes().at(0);
	std::set < Node_ptr > iNodeSet; //!< InitialNodeSet for DFA

	for (std::vector<Node_ptr>::size_type i = 0; i < NodeList.size(); ++i) {

		std::queue < Node_ptr > NodeQueue;

		/** @note temp is the computed Null-closure for NodeList[i] */
		std::set < Node_ptr > temp;

		NodeQueue.push(NodeList[i]);

		while (!NodeQueue.empty()) {

			std::set < Node_ptr > scanned_nodes;

			Node_ptr unaryNode = NodeQueue.front();

			scanned_nodes.insert(unaryNode);
			NodeQueue.pop();

			globalmapIterator = globalmap.find(unaryNode);

			if (globalmapIterator != globalmap.end())
				mymapIterator = globalmapIterator->second.find(
						label_cache::label_id_type(0));

			if (mymapIterator != globalmapIterator->second.end()) {
				std::set < Node_ptr > temp_set = (*mymapIterator).second;

				for (std::set<Node_ptr>::iterator it = temp_set.begin(); it
						!= temp_set.end(); it++) {
					Node_ptr unaryNode;
					unaryNode = *it;
					if (!(scanned_nodes.count(unaryNode) > 0)) {
						temp.insert(unaryNode);
						NodeQueue.push(unaryNode);
						scanned_nodes.insert(unaryNode);
					}
				} //!< for loop ends
			} //!< if(mymapIterator != globalmapIterator->second.end())
			else {
				//!< do nothing since there is no null transition from unaryNode
			}

		} //!< while(!NodeQueue.empty())

		globalmapIterator = globalmap.find(NodeList[i]);

		mymapIterator = globalmapIterator->second.find(
				label_cache::label_id_type(0));

		if (!temp.empty() && mymapIterator != globalmapIterator->second.end()) {
			(*mymapIterator).second.insert(temp.begin(), temp.end());
		}

		/** If current Node is the initial Node, then modify initial NodeSet */
		if (iNode == NodeList[i]) {
			iNodeSet.insert(NodeList[i]);
			if (mymapIterator != globalmapIterator->second.end()) {
				std::set < Node_ptr > mymapSet = (*mymapIterator).second;
				iNodeSet.insert(mymapSet.begin(), mymapSet.end());
			}
		}

	} //!< Null closures computed

	/** Merge Null closures into the sets. */

	for (std::vector<Node_ptr>::size_type i = 0; i < NodeList.size(); ++i) {

		std::vector < AutomatonEdge > myEdges;

		myEdges = NodeList[i]->getEdges();

		globalmapIterator = globalmap.find(NodeList[i]);

		if (globalmapIterator == globalmap.end())
			throw std::runtime_error("The node is not there in the globalmap. ");

		for (std::vector<AutomatonEdge>::size_type j = 0; j < myEdges.size(); ++j) {

			label_cache::label_id_type lId = myEdges[j].getEdgeLabelId();

			if (lId != label_cache::label_id_type(0)) {

				std::set < Node_ptr > nullClosure_set;

				mymapIterator = globalmapIterator->second.find(lId);

				std::set < Node_ptr > aNodeSet = mymapIterator->second;

				std::map<Node_ptr, std::map<label_cache::label_id_type,
						std::set<Node_ptr> > >::iterator temp_globalmapIterator;
				for (std::set<Node_ptr>::iterator nIt = aNodeSet.begin(); nIt
						!= aNodeSet.end(); ++nIt) {
					temp_globalmapIterator = globalmap.find(*nIt);
					std::map < label_cache::label_id_type, std::set<Node_ptr>
							> temp_mymap;
					std::map<label_cache::label_id_type, std::set<Node_ptr> >::iterator
							temp_mymapIterator;
					temp_mymap = temp_globalmapIterator->second;
					temp_mymapIterator = temp_mymap.find(
							label_cache::label_id_type(0));
					if (temp_mymapIterator != temp_mymap.end()) {
						std::set < Node_ptr > temp_nullClosure_set
								= temp_mymapIterator->second;
						nullClosure_set.insert(temp_nullClosure_set.begin(),
								temp_nullClosure_set.end());
					}
				}
				if (!nullClosure_set.empty())
					(*mymapIterator).second.insert(nullClosure_set.begin(),
							nullClosure_set.end());
			}
		}
	}

	/** Create a set from all NodeList */
	std::set < Node_ptr > tempNodeSet;
	for (std::vector<Node_ptr>::size_type i = 0; i < NodeList.size(); ++i)
		tempNodeSet.insert(NodeList[i]);

	Node_ptr* NodeArray;
	int NumberofNodes = tempNodeSet.size();
	NodeArray = new Node_ptr[NumberofNodes];

	/** Create a NodeArray of type T from tempNodeSet */
	int i = 0;
	for (std::set<Node_ptr>::iterator it = tempNodeSet.begin(); it
			!= tempNodeSet.end(); ++it, ++i)
		NodeArray[i] = *it;

	/**
	 * Create mySet from NodeArray.
	 *
	 * @note mySet is used to compute subsets of a power-set
	 */
	Theset < Node_ptr > mySet(NodeArray, NumberofNodes);

	/** Maps a newly constructed subset to its corresponding node
	 *
	 * @note Used later while adding adjacent nodes in a DFA
	 */
	std::map < std::set<Node_ptr>, Node_ptr > setToNodemap;

	DFAutomaton DFAOfNFA; //!< Create DFA

	/** Vector contains all subsets of NFA nodes' ids. */
	std::vector < std::set<Node_ptr> > DFANodesetV;

	NodeList.clear();

	/** Compute Power set. */
	for (int i = 0; i < (1 << NumberofNodes); i++) {

		std::set < Node_ptr > NewSubset = mySet.NextSubset();
		const AutomatonNode& newNodeObj_ref = DFAOfNFA.createNode();
		Node_ptr newNode = access_method::access(newNodeObj_ref);
		NodeList.push_back(newNode);
		DFANodesetV.push_back(NewSubset);
	}

	Node_ptr DFA_initialNode;
	for (std::vector<std::set<Node_ptr> >::size_type i = 0; i
			< DFANodesetV.size(); ++i) {

		/** There is only one initial state in a DFA */
		if (DFANodesetV[i] == iNodeSet)
			DFA_initialNode = NodeList[i];
	}
	DFAOfNFA.setInitialNode(*DFA_initialNode);

	std::vector < Node_ptr > NFAacceptNodes = thisNFA_copy.getAcceptNodes();

	/** Compute DFA transitions as well as reachable set. */
	for (std::vector<Node_ptr>::size_type rSize = 0; rSize < NodeList.size(); rSize++) {

		Node_ptr rNode = NodeList.at(rSize);
		std::set < Node_ptr > rNodeSet = DFANodesetV.at(rSize);

		/** Find accept nodes of DFA */
		bool accept_node_found = false;

		/** for each NFAacceptNodes */
		for (std::vector<Node_ptr>::iterator Nit = NFAacceptNodes.begin(); !accept_node_found
				&& Nit != NFAacceptNodes.end(); ++Nit) {

			Node_ptr NFA_acceptNode = *Nit;

			/** for each element of rNodeSet */
			for (std::set<Node_ptr>::iterator it = rNodeSet.begin(); it
					!= rNodeSet.end(); it++) {

				/** exit the loop as soon as an accepting nId is found in rNodeSet */
				if (NFA_acceptNode == *it) {
					accept_node_found = true;
					DFAOfNFA.setAcceptNode(*rNode);
					break;
				}
			}
		} //!< end-Find accept nodes of DFA

		std::set < label_cache::label_id_type > labelIds
				= thisNFA_copy.getAlphabetIds();

		/** for each label */
		for (std::set<label_cache::label_id_type>::iterator lIt =
				labelIds.begin(); lIt != labelIds.end(); ++lIt) {

			/** Only for non-EPS labels */
			if (*lIt != label_cache::label_id_type(0)) {

				/** adjacentNodeSet for each DFA node */
				std::set < Node_ptr > adjacentNodeSet;

				/** Compute the adjacent node from NFA for each element node of a DFA node
				 * and add it to the adjacentNodeSet. */
				for (std::set<Node_ptr>::iterator it = rNodeSet.begin(); it
						!= rNodeSet.end(); it++) {
					Node_ptr unaryNode = *it;
					globalmapIterator = globalmap.find(unaryNode);
					if (globalmapIterator != globalmap.end()) {

						/** Check for label k */
						mymapIterator
								= ((*globalmapIterator).second.find(*lIt));

						if (mymapIterator != (*globalmapIterator).second.end()) {
							adjacentNodeSet.insert(
									(*mymapIterator).second.begin(),
									(*mymapIterator).second.end());

						} //!< if(mymapIterator != (*globalmapIterator).second.end())
					} //!< if(globalmapIterator != globalmap.end())
				} //!< adjacent Node computed

				/** find the Index of computed adjacentNodeset in DFANodesetV */
				std::vector<std::set<Node_ptr> >::size_type adjNodeIdx = 0;
				for (std::vector<std::set<Node_ptr> >::iterator it =
						DFANodesetV.begin(); it != DFANodesetV.end(); ++it, ++adjNodeIdx)
					if (*it == adjacentNodeSet) {
						break;
					}

				/** Check if the adjacentNode is part of DFANodeset Vector */
				if (adjNodeIdx < DFANodesetV.size()) {

					std::string myLabel = label_cache::getLabel(*lIt);

					DFAOfNFA.addToAlphabet(myLabel);
					rNode->addSuccNode(*(NodeList.at(adjNodeIdx)), 1, *lIt);

				} //!< if(adjNodeIdx < DFANodesetV.size())
				else
					throw std::runtime_error(
							"Adjacent Node computed is not part of DFANodeset Vector");
			} //!< if(label_Id[k] != 0)
		} //!< for each label

	} //!< while(!reachable_nodes_queue.empty())

	/** Compute reachable set of Nodes. */
	std::vector < Node_ptr > reachDFANodeList;
	std::queue < Node_ptr > reachDFANodeQueue;
	DFA_initialNode->setStatus(AutomatonNode::IS_VISITED);
	reachDFANodeList.push_back(DFA_initialNode);
	reachDFANodeQueue.push(DFA_initialNode);
	while (!reachDFANodeQueue.empty()) {
		Node_ptr rNode = reachDFANodeQueue.front();
		reachDFANodeQueue.pop();
		std::vector < AutomatonEdge > myEdges = rNode->getEdges();
		for (std::vector<AutomatonEdge>::iterator eIt = myEdges.begin(); eIt
				!= myEdges.end(); ++eIt) {
			const AutomatonNode& succNode_ref = (*eIt).getSuccessorNode();
			Node_ptr succNode = access_method::access(succNode_ref);
			if (succNode->getStatus() == AutomatonNode::IS_NOT_VISITED) {
				succNode->setStatus(AutomatonNode::IS_VISITED);
				reachDFANodeList.push_back(succNode);
				reachDFANodeQueue.push(succNode);
			}
		}
	}

	/** Reset nodes' statuses. */
	for (std::vector<Node_ptr>::iterator nIt = reachDFANodeList.begin(); nIt
			!= reachDFANodeList.end(); ++nIt)
		(*nIt)->setStatus(AutomatonNode::IS_NOT_VISITED);

	std::vector < new_node_id::node_id_type > NodeIds = DFAOfNFA.getNodeIds();
	NodeIds.resize(reachDFANodeList.size());

	DFAOfNFA.resetNodeList(reachDFANodeList);
	DFAOfNFA.resetNodeIds(NodeIds);

	// minimize DFA
	DFAutomaton minDFA = DFAOfNFA.minimize();

	if (!(minDFA.is_empty()))
		return minDFA;
	else
		throw std::runtime_error("Minimized DFA is empty.");
}

std::vector<std::string> NFAutomaton::getShortestPath() const {

	NFAutomaton thisAutomaton_copy(*this);
	std::vector < Node_ptr > myInitialNodes
			= thisAutomaton_copy.getInitialNodes();

	if (myInitialNodes.size() > 1)
		thisAutomaton_copy.mergeInitialNodes(false);

	Node_ptr source = thisAutomaton_copy.getInitialNodes().at(0);
	std::vector < Node_ptr > myAcceptNodes
			= thisAutomaton_copy.getAcceptNodes();
	distance_type max_distance = distance_type(MAX_DISTANCE);

	std::vector < Node_ptr > NodeList = thisAutomaton_copy.getNodeList();

	std::map < Node_ptr, temp_info_type > temp_info_map;
	std::map<Node_ptr, temp_info_type>::iterator temp_info_mapIt;

	for (std::vector<Node_ptr>::iterator nIt = NodeList.begin(); nIt
			!= NodeList.end(); ++nIt) {

		temp_info_type temp_info;

		if (*nIt == source)
			temp_info.distance = distance_type(0);
		else
			temp_info.distance = max_distance;

		temp_info_map.insert(std::pair<Node_ptr, temp_info_type>(*nIt,
				temp_info));
	}

	std::queue < Node_ptr > reachable_nodes_queue;
	reachable_nodes_queue.push(source);

	while (!reachable_nodes_queue.empty()) {

		Node_ptr rNode = reachable_nodes_queue.front();
		reachable_nodes_queue.pop();

		distance_type u_distance = distance_type(0);
		temp_info_mapIt = temp_info_map.find(rNode);

		if (temp_info_mapIt != temp_info_map.end()) {
			temp_info_type temp_info = temp_info_mapIt->second;
			u_distance = temp_info.distance;
		}

		std::vector < AutomatonEdge > myEdges = rNode->getEdges();

		for (std::vector<AutomatonEdge>::iterator eIt = myEdges.begin(); eIt
				!= myEdges.end(); ++eIt) {

			distance_type v_distance = distance_type(0);

			distance_type uv_weight = (*eIt).getEdgeCost();
			const AutomatonNode& succNodeObj_ref = (*eIt).getSuccessorNode();
			Node_ptr mySuccNode = access_method::access(succNodeObj_ref);
			label_cache::label_id_type lId = (*eIt).getEdgeLabelId();

			temp_info_mapIt = temp_info_map.find(mySuccNode);
			if (temp_info_mapIt != temp_info_map.end()) {
				temp_info_type temp_info = temp_info_mapIt->second;
				v_distance = temp_info.distance;
			}

			if ((u_distance + uv_weight) < v_distance) {
				(temp_info_mapIt->second).distance = (u_distance + uv_weight);
				(temp_info_mapIt->second).parent = rNode;
				(temp_info_mapIt->second).label_id = lId;
				reachable_nodes_queue.push(mySuccNode);
			}
		}
	}

	distance_type shortest_distance = max_distance;
	Node_ptr shortest_accept_node;

	/** Get the distance for an accepting state with least distance. */
	for (std::vector<Node_ptr>::iterator aIt = myAcceptNodes.begin(); aIt
			!= myAcceptNodes.end(); ++aIt) {
		distance_type temp_distance = max_distance;
		temp_info_mapIt = temp_info_map.find(*aIt);
		if (temp_info_mapIt != temp_info_map.end()) {
			temp_info_type temp_info = temp_info_mapIt->second;
			temp_distance = temp_info.distance;
		}
		if (temp_distance < shortest_distance) {
			shortest_distance = temp_distance;
			shortest_accept_node = *aIt;
		}
	}

	/** Revert back to the original configuration. */
	if (myInitialNodes.size() > 1)
		thisAutomaton_copy.separateInitialNodes();

	if (shortest_distance == max_distance) {
		throw std::runtime_error(
				"@ShortestPath: There is no accepting state in the automaton.");
	}

	std::vector < std::string > shortest_path;

	/** Get the path labels */
	if (shortest_accept_node != source) {

		Node_ptr parent = shortest_accept_node;
		std::vector < label_cache::label_id_type > Rev_shortest_path;
		temp_info_type temp_info;

		temp_info_mapIt = temp_info_map.find(parent);
		if (temp_info_mapIt != temp_info_map.end()) {
			temp_info = temp_info_mapIt->second;
			parent = temp_info.parent;
			label_cache::label_id_type label_id = temp_info.label_id;
			Rev_shortest_path.push_back(label_id);
		}

		while (parent != source) {
			temp_info_mapIt = temp_info_map.find(parent);
			if (temp_info_mapIt != temp_info_map.end()) {
				temp_info = temp_info_mapIt->second;
				parent = temp_info.parent;
				label_cache::label_id_type label_id = temp_info.label_id;
				Rev_shortest_path.push_back(label_id);
			}
		}

		if (myInitialNodes.size() > 1)
			Rev_shortest_path.erase(Rev_shortest_path.begin()
					+ (Rev_shortest_path.size() - 1));

		std::vector<label_cache::label_id_type>::reverse_iterator rIt;
		for (rIt = Rev_shortest_path.rbegin(); rIt < Rev_shortest_path.rend(); ++rIt) {
			std::string label = label_cache::getLabel(*rIt);
			shortest_path.push_back(label);
		}

	}
	return shortest_path;
}
}
#endif /* NFAUTOMATON_HPP_ */
