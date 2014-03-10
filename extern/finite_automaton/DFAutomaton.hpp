#ifndef DFAUTOMATON_HPP_
#define DFAUTOMATON_HPP_

#include "DFAutomaton.h"

namespace finite_automaton {

/**
 * Given an alphabet (set of labels), the method creates
 * a chaos.
 *
 * Chaos is a language that accepts all the possible words
 * created of the labels from its alphabet. */
DFAutomaton createChaos(const std::set<std::string>& alphabet) {

	DFAutomaton Chaos;

	const AutomatonNode& newNode = Chaos.createNode();
	Chaos.setAcceptNode(newNode);
	Chaos.setInitialNode(newNode);

	for (std::set<std::string>::iterator it = alphabet.begin(); it
			!= alphabet.end(); ++it) {
		Chaos.addToAlphabet(*it);
		Chaos.addEdge(newNode, newNode, 1, *it);
	}

	return Chaos;
}

/** Given an alphabet (set of labels), the method creates
 * an empty language.
 *
 * It is an empty DFA with its alphabet set. */
DFAutomaton createEmptyLanguage(const std::set<std::string>& alphabet) {
	DFAutomaton EmptyLang;
	EmptyLang.setAlphabet(alphabet);
	return EmptyLang;
}

/** Given a string, the method creates a DFA which determines
 * a language recognizing the given string. */
DFAutomaton createDFAFromString(const std::string& str) {

	DFAutomaton strDFA;

	const AutomatonNode& newNode = strDFA.createNode();
	strDFA.setInitialNode(newNode);
	Node_ptr prevNode = access_method::access(newNode);

	for (unsigned int i = 0; i < str.size(); ++i) {
		const AutomatonNode& nextNode = strDFA.createNode();
		std::string mystring;
		mystring = mystring + str.at(i);
		strDFA.addToAlphabet(mystring);
		strDFA.addEdge(*prevNode, nextNode, 1, mystring);
		prevNode = access_method::access(nextNode);
	}
	strDFA.setAcceptNode(*prevNode);
	return strDFA;
}

/** Given a vector of labels i.e., a word, the method creates
 * a DFA. This DFA determines a language that accepts this word. */
DFAutomaton createDFAFromWord(const std::vector<std::string>& word) {

	DFAutomaton wordDFA;

	AutomatonNode newNode = wordDFA.createNode();
	wordDFA.setInitialNode(newNode);
	AutomatonNode prevNode = newNode;

	for (std::vector<std::string>::const_iterator sIt = word.begin(); sIt
			!= word.end(); ++sIt) {

		AutomatonNode nextNode = wordDFA.createNode();
		wordDFA.addToAlphabet(*sIt);
		wordDFA.addEdge(prevNode, nextNode, 1, *sIt);
		prevNode = nextNode;
	}
	wordDFA.setAcceptNode(prevNode);
	return wordDFA;
}

DFAutomaton::DFAutomaton() {
}

DFAutomaton::~DFAutomaton() {
}

DFAutomaton::DFAutomaton(const DFAutomaton& automaton) {

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

NFAutomaton DFAutomaton::convertToNFAStar() const {

	NFAutomaton NFAStar = this->convertToNFA();
	NFAStar.addToAlphabet("EPS");
	std::vector < Node_ptr > NFAStar_NodeList = NFAStar.getNodeList();

	for (std::vector<Node_ptr>::iterator nIt = NFAStar_NodeList.begin(); nIt
			!= NFAStar_NodeList.end(); ++nIt)
		NFAStar.addEdge(*(*nIt), *(*nIt), 1, "EPS");

	return NFAStar;
}

bool DFAutomaton::Accepts(const std::string word) const {

	std::vector < Node_ptr > myAcceptNodes = this->getAcceptNodes();

	Node_ptr currentNode = this->getInitialNodes().at(0);

	for (unsigned int i = 0; i < word.size(); ++i) {
		std::string mystring;
		mystring = mystring + word.at(i);

		std::vector < AutomatonEdge > myEdges = currentNode->getEdges();

		for (std::vector<AutomatonEdge>::iterator eIt = myEdges.begin(); eIt
				!= myEdges.end(); ++eIt) {

			std::string myLabel =
					label_cache::getLabel((*eIt).getEdgeLabelId());
			if (myLabel.compare(mystring) == 0) {
				const AutomatonNode& succNodeObj_ref =
						(*eIt).getSuccessorNode();
				currentNode = access_method::access(succNodeObj_ref);
				break;
			}
		}
	}

	for (std::vector<Node_ptr>::iterator aIt = myAcceptNodes.begin(); aIt
			!= myAcceptNodes.end(); ++aIt) {
		if (currentNode == (*aIt))
			return true;
	}

	return false;
}

DFAutomaton DFAutomaton::concatenate(const DFAutomaton& secondDFA) const {

	DFAutomaton firstDFA_copy(*this);
	DFAutomaton secondDFA_copy(secondDFA);

	std::set < label_cache::label_id_type > firstDFA_labelIds
			= firstDFA_copy.getAlphabetIds();
	std::set < label_cache::label_id_type > secondDFA_labelIds
			= secondDFA_copy.getAlphabetIds();

	firstDFA_copy.extend(secondDFA_labelIds);
	secondDFA_copy.extend(firstDFA_labelIds);

	NFAutomaton firstNFA = firstDFA_copy.convertToNFA();
	NFAutomaton secondNFA = secondDFA_copy.convertToNFA();
	NFAutomaton concatenated_NFA = firstNFA.concatenate(secondNFA);
	DFAutomaton concatenated_DFA = concatenated_NFA.convertToDFA();
	return concatenated_DFA;
}

NFAutomaton DFAutomaton::convertToNFA() const {

	NFAutomaton NFAOfDFA;
	std::vector < Node_ptr > myAcceptNodes = getAcceptNodes();
	std::set < Node_ptr > myAcceptNodesSet;
	for (std::vector<Node_ptr>::iterator aIt = myAcceptNodes.begin(); aIt
			!= myAcceptNodes.end(); ++aIt)
		myAcceptNodesSet.insert(*aIt);

	std::vector < Node_ptr > myInitialNodes = getInitialNodes();

	std::vector < Node_ptr > myNodes = getNodeList();

	std::map < Node_ptr, Node_ptr > origNodetoNFANodemap;
	std::map<Node_ptr, Node_ptr>::iterator origNodetoNFANodemapIt;

	std::vector < Node_ptr > NFANodes;

	/** Copy nodes including initial and accept nodes */
	for (std::vector<Node_ptr>::size_type i = 0; i < myNodes.size(); ++i) {
		const AutomatonNode& newNodeObj_ref = NFAOfDFA.createNode();
		Node_ptr newNode = access_method::access(newNodeObj_ref);
		NFANodes.push_back(newNode);
		if (myAcceptNodesSet.count(myNodes[i]) > 0)
			NFAOfDFA.setAcceptNode(*newNode);
		if (myInitialNodes[0] == myNodes[i])
			NFAOfDFA.setInitialNode(*newNode);
		origNodetoNFANodemap.insert(std::pair<Node_ptr, Node_ptr>(myNodes[i],
				newNode));
	}

	/** Copy transitions */
	for (std::vector<Node_ptr>::size_type i = 0; i < myNodes.size(); ++i) {
		std::vector < AutomatonEdge > myEdges = myNodes[i]->getEdges();

		for (std::vector<AutomatonEdge>::size_type j = 0; j < myEdges.size(); ++j) {
			const AutomatonNode& succNodeObj_ref =
					myEdges[j].getSuccessorNode();
			origNodetoNFANodemapIt = origNodetoNFANodemap.find(
					access_method::access(succNodeObj_ref));

			if (origNodetoNFANodemapIt != origNodetoNFANodemap.end()) {

				std::string myLabel = label_cache::getLabel(
						myEdges[j].getEdgeLabelId());
				NFAOfDFA.addToAlphabet(myLabel);
				NFANodes[i]->addSuccNode(*(origNodetoNFANodemapIt->second),
						myEdges[j].getEdgeCost(), myEdges[j].getEdgeLabelId());
			}
		}
	}
	return NFAOfDFA;
}

void DFAutomaton::makeComplete() {

	/** label id's in an alphabet are mapped to their corresponding
	 * indices in the alphabet. */
	std::set < label_cache::label_id_type > labelIds = this->getAlphabetIds();
	std::map<label_cache::label_id_type, int> label_id_to_index_map;
	std::map<label_cache::label_id_type, int>::iterator label_id_to_index_mapIt;
	std::set<label_cache::label_id_type>::size_type i = 0;
	for (std::set<label_cache::label_id_type>::iterator lIt = labelIds.begin(); lIt
			!= labelIds.end(); ++lIt, ++i) {
		label_id_to_index_map.insert(
				std::pair<label_cache::label_id_type, int>(*lIt, i));
	}

	/** Create a dead node. */
	Node_ptr deadNode;
	bool node_added = false;

	std::vector < Node_ptr > NodeList = this->getNodeList();

	for (std::vector<Node_ptr>::iterator nIt = NodeList.begin(); nIt
			!= NodeList.end(); ++nIt) {

		std::vector < AutomatonEdge > myEdges = (*nIt)->getEdges();

		/** Only if an edge is missing for some label */
		if (myEdges.size() != labelIds.size()) {

			/** dead node is added to the graph only once. */
			if (!node_added) {

				/** dead node is added to the automaton. */
				const AutomatonNode& newNodeObj_ref = this->createNode();
				deadNode = access_method::access(newNodeObj_ref);

				/** Add transitions from dead node to itself for each alphabet label. */
				for (std::set<label_cache::label_id_type>::iterator lIt =
						labelIds.begin(); lIt != labelIds.end(); ++lIt)
					deadNode->addSuccNode(*deadNode, 1, *lIt);

				node_added = true; //!< set the flag
			}

			/** edge_available vector helps to determine if each DFA node
			 * has an edge for every label. */
			std::vector<bool> edge_available(labelIds.size());

			/** scan all edges of a node and make entries in edge_available for those edges. */
			for (std::vector<AutomatonEdge>::size_type j = 0; j
					< myEdges.size(); ++j) {
				label_cache::label_id_type label =
						myEdges.at(j).getEdgeLabelId();
				label_id_to_index_mapIt = label_id_to_index_map.find(label);
				if (label_id_to_index_mapIt != label_id_to_index_map.end()) {
					edge_available.at(label_id_to_index_mapIt->second) = true;
				}
			}

			/** Add edges for all missing labels for this node. */

			std::set<label_cache::label_id_type>::size_type k = 0;
			for (std::set<label_cache::label_id_type>::iterator lIt =
					labelIds.begin(); lIt != labelIds.end(); ++lIt, ++k) {

				/** only for the absent label edges. */
				if (edge_available.at(k) == false) {
					(*nIt)->addSuccNode(*deadNode, 1, *lIt);
				}
			} //!< end-for labelIds
		} //!< end if(myEdges.size() != labelIds.size())
	} //!< end-for Nodelist

	if (!node_added) {
		//!< The graph is already complete;
	}
	return;
}

void DFAutomaton::extend(const std::set<label_cache::label_id_type> LabelIds) {

	std::set < label_cache::label_id_type > labelIdSet = this->getAlphabetIds();

	for (std::set<label_cache::label_id_type>::iterator lIt = LabelIds.begin(); lIt
			!= LabelIds.end(); ++lIt)
		labelIdSet.insert(*lIt);

	this->setAlphabetIds(labelIdSet);
	this->makeComplete();

	return;
}

NFAutomaton DFAutomaton::computeKleeneStar() const {

	std::vector < Node_ptr > myAcceptNodes = getAcceptNodes();
	std::set < Node_ptr > myAcceptNodesSet;
	for (std::vector<Node_ptr>::iterator aIt = myAcceptNodes.begin(); aIt
			!= myAcceptNodes.end(); ++aIt)
		myAcceptNodesSet.insert(*aIt);

	std::vector < Node_ptr > myInitialNodes = getInitialNodes();

	std::vector < Node_ptr > myNodes = getNodeList();

	std::map < Node_ptr, Node_ptr > origNodetoNFANodemap;
	std::map<Node_ptr, Node_ptr>::iterator origNodetoNFANodemapIt;

	NFAutomaton NFAOfDFA;
	std::vector < Node_ptr > NFANodes;
	Node_ptr oldInitialNode;
	myAcceptNodes.clear();

	/** Copy nodes including initial and accept nodes */
	for (std::vector<Node_ptr>::size_type i = 0; i < myNodes.size(); ++i) {
		const AutomatonNode& newNodeObj_ref = NFAOfDFA.createNode();
		Node_ptr newNode = access_method::access(newNodeObj_ref);
		NFANodes.push_back(newNode);
		if (myAcceptNodesSet.count(myNodes[i]) > 0) {
			myAcceptNodes.push_back(newNode);
			NFAOfDFA.setAcceptNode(*newNode);
		}
		if (myInitialNodes[0] == myNodes[i])
			oldInitialNode = newNode;
		origNodetoNFANodemap.insert(std::pair<Node_ptr, Node_ptr>(myNodes[i],
				newNode));
	}

	/** Copy transitions */
	for (std::vector<Node_ptr>::size_type i = 0; i < myNodes.size(); ++i) {
		std::vector < AutomatonEdge > myEdges = myNodes[i]->getEdges();

		for (std::vector<AutomatonEdge>::size_type j = 0; j < myEdges.size(); ++j) {
			const AutomatonNode& succNodeObj_ref =
					myEdges[j].getSuccessorNode();
			origNodetoNFANodemapIt = origNodetoNFANodemap.find(
					access_method::access(succNodeObj_ref));

			if (origNodetoNFANodemapIt != origNodetoNFANodemap.end()) {

				std::string myLabel = label_cache::getLabel(
						myEdges[j].getEdgeLabelId());
				NFAOfDFA.addToAlphabet(myLabel);
				NFANodes[i]->addSuccNode(*(origNodetoNFANodemapIt->second),
						myEdges[j].getEdgeCost(), myEdges[j].getEdgeLabelId());
			}
		}
	}

	NFAOfDFA.addToAlphabet("EPS");

	/** Add Null transitions from accept states to oldInitialState */
	for (std::vector<Node_ptr>::size_type i = 0; i < myAcceptNodes.size(); ++i)
		myAcceptNodes[i]->addSuccNode(*oldInitialNode, 1,
				label_cache::label_id_type(0));

	const AutomatonNode& newNodeObj_ref = NFAOfDFA.createNode();
	Node_ptr newNode = access_method::access(newNodeObj_ref);

	newNode->addSuccNode(*oldInitialNode, 1, label_cache::label_id_type(0));

	/** make newNode initial node */
	NFAOfDFA.setInitialNode(*newNode);

	/** make newNode accepting node */
	NFAOfDFA.setAcceptNode(*newNode);

	return NFAOfDFA;
}

DFAutomaton DFAutomaton::computeDifference(const DFAutomaton& secondDFA) const {

	DFAutomaton firstDFA_copy(*this);
	DFAutomaton secondDFA_copy(secondDFA);

	std::set < label_cache::label_id_type > second_labelIdSet
			= secondDFA_copy.getAlphabetIds();
	std::set < label_cache::label_id_type > first_labelIdSet
			= firstDFA_copy.getAlphabetIds();

	secondDFA_copy.extend(first_labelIdSet);
	firstDFA_copy.extend(second_labelIdSet);

	first_labelIdSet.clear();
	second_labelIdSet.clear();

	/** If minuend is empty, then difference is also empty */
	if (this->getNodeList().empty())
		return *this;

	DFAutomaton complementOfsecondDFA = secondDFA_copy.computeComplement();

	DFAutomaton DifferenceAutomaton = firstDFA_copy.computeIntersection(
			complementOfsecondDFA);
	return DifferenceAutomaton;
}

DFAutomaton DFAutomaton::computeComplement() const {

	/** The complement of an empty language is a Chaos. */
	if (getNodeList().empty()) {
		std::set < std::string > alphabet = getAlphabet();
		DFAutomaton this_chaos = createChaos(alphabet);
		return this_chaos;
	}

	DFAutomaton thisDFA_copy(*this);

	std::map < Node_ptr, Node_ptr > origNodetoCompNodemap;
	std::map<Node_ptr, Node_ptr>::iterator origNodetoCompNodemapIt;

	std::vector < Node_ptr > myAcceptNodesList = thisDFA_copy.getAcceptNodes();
	std::set < Node_ptr > myAcceptNodesSet;
	for (std::vector<Node_ptr>::iterator aIt = myAcceptNodesList.begin(); aIt
			!= myAcceptNodesList.end(); ++aIt)
		myAcceptNodesSet.insert(*aIt);

	std::vector < Node_ptr > myInitialNode = thisDFA_copy.getInitialNodes();
	std::vector < Node_ptr > myNodeList = thisDFA_copy.getNodeList();

	DFAutomaton myComplement;

	for (std::vector<Node_ptr>::size_type i = 0; i < myNodeList.size(); ++i) {
		const AutomatonNode& newNodeObj_ref = myComplement.createNode();
		Node_ptr newNode = access_method::access(newNodeObj_ref);
		if (!(myAcceptNodesSet.count(myNodeList[i]) > 0))
			myComplement.setAcceptNode(*newNode);
		if (myNodeList[i] == myInitialNode.at(0))
			myComplement.setInitialNode(*newNode);
		origNodetoCompNodemap.insert(std::pair<Node_ptr, Node_ptr>(
				myNodeList[i], newNode));
	}

	std::set < label_cache::label_id_type > labelIds
			= thisDFA_copy.getAlphabetIds();
	myComplement.setAlphabetIds(labelIds);

	std::vector < Node_ptr > CompNodes = myComplement.getNodeList();

	/** Copy transitions */
	for (std::vector<Node_ptr>::size_type i = 0; i < myNodeList.size(); ++i) {
		std::vector < AutomatonEdge > myEdges = myNodeList[i]->getEdges();

		for (std::vector<AutomatonEdge>::size_type j = 0; j < myEdges.size(); ++j) {
			const AutomatonNode& succNodeObj_ref =
					myEdges[j].getSuccessorNode();
			origNodetoCompNodemapIt = origNodetoCompNodemap.find(
					access_method::access(succNodeObj_ref));

			if (origNodetoCompNodemapIt != origNodetoCompNodemap.end())
				CompNodes[i]->addSuccNode(*(origNodetoCompNodemapIt->second),
						myEdges[j].getEdgeCost(), myEdges[j].getEdgeLabelId());
		}
	}
	return myComplement;
}

DFAutomaton DFAutomaton::computeUnion(const DFAutomaton& secondDFA) const {

	DFAutomaton firstDFA_copy(*this);
	DFAutomaton secondDFA_copy(secondDFA);

	std::set < label_cache::label_id_type > second_labelIdSet
			= secondDFA_copy.getAlphabetIds();
	std::set < label_cache::label_id_type > first_labelIdSet
			= firstDFA_copy.getAlphabetIds();

	firstDFA_copy.extend(second_labelIdSet);
	secondDFA_copy.extend(first_labelIdSet);

	first_labelIdSet.clear();
	second_labelIdSet.clear();

	/** If any of the DFA is empty, return the other one as Union. */
	if (firstDFA_copy.getNodeList().empty())
		return secondDFA_copy;
	else if (secondDFA.getNodeList().empty())
		return firstDFA_copy;

	/** label ids in an alphabet are mapped to their corresponding
	 * indices in the alphabet. */
	std::set < label_cache::label_id_type > labelIds
			= firstDFA_copy.getAlphabetIds();
	std::map<label_cache::label_id_type, unsigned int> label_id_to_index_map;
	std::map<label_cache::label_id_type, unsigned int>::iterator
			label_id_to_index_map_It;
	std::set<label_cache::label_id_type>::size_type i = 0;
	for (std::set<label_cache::label_id_type>::iterator lIt = labelIds.begin(); lIt
			!= labelIds.end(); ++lIt, ++i) {
		label_id_to_index_map.insert(std::pair<label_cache::label_id_type,
				unsigned int>(*lIt, i));
	}

	/** DFAmap - maps each DFA node to a vector of nodes.
	 *
	 *  For a node, entry at Index i in its corresponding vector
	 *  represents the successor node reached from transition on label
	 *  at index i in DFAalphabet.
	 */
	std::map<Node_ptr, std::vector<Node_ptr> > thisDFAmap, secondDFAmap;
	std::map<Node_ptr, std::vector<Node_ptr> >::iterator DFAmapIterator1,
			DFAmapIterator2;

	std::vector < Node_ptr > thisNodeList = firstDFA_copy.getNodeList();

	for (std::vector<Node_ptr>::iterator it = thisNodeList.begin(); it
			!= thisNodeList.end(); ++it) {

		/** Initialize DFAmap for each node with an emptyDFANodeVector
		 * to be filled later. */
		std::vector < Node_ptr > emptyDFANodeVector(labelIds.size());
		thisDFAmap.insert(std::pair<Node_ptr, std::vector<Node_ptr> >(*it,
				emptyDFANodeVector));

		DFAmapIterator1 = thisDFAmap.find(*it);

		std::vector < AutomatonEdge > myEdges = (*it)->getEdges();

		/** fill emptyDFANodeVector */
		for (std::vector<AutomatonEdge>::iterator edgeIt = myEdges.begin(); edgeIt
				!= myEdges.end(); ++edgeIt) {

			const AutomatonNode& succNodeObj_ref = edgeIt->getSuccessorNode();
			Node_ptr succNode = access_method::access(succNodeObj_ref);
			label_id_to_index_map_It = label_id_to_index_map.find(
					edgeIt->getEdgeLabelId());
			if (label_id_to_index_map_It != label_id_to_index_map.end())
				DFAmapIterator1->second.at(label_id_to_index_map_It->second)
						= succNode;
		}
	}

	/** get initial node and accept nodes */
	Node_ptr thisinitialNode = firstDFA_copy.getInitialNodes().at(0);
	std::set < Node_ptr > thisacceptNodesSet;
	std::vector < Node_ptr > thisacceptNodeList
			= firstDFA_copy.getAcceptNodes();
	for (std::vector<Node_ptr>::iterator it = thisacceptNodeList.begin(); it
			!= thisacceptNodeList.end(); ++it)
		thisacceptNodesSet.insert(*it);

	/** Repeat the same for secondDFA */
	std::vector < Node_ptr > secondDFANodeList = secondDFA_copy.getNodeList();

	for (std::vector<Node_ptr>::iterator it = secondDFANodeList.begin(); it
			!= secondDFANodeList.end(); ++it) {

		/** Initialize DFAmap for each node with an emptyDFANodeVector, which
		 * will be filled later. */
		std::vector < Node_ptr > emptyDFANodeVector(labelIds.size());
		secondDFAmap.insert(std::pair<Node_ptr, std::vector<Node_ptr> >(*it,
				emptyDFANodeVector));

		DFAmapIterator2 = secondDFAmap.find(*it);

		if (DFAmapIterator2 == secondDFAmap.end())
			throw std::runtime_error(
					"@Union: Node is not there in secondDFAmap.");

		std::vector < AutomatonEdge > myEdges = (*it)->getEdges();

		for (std::vector<AutomatonEdge>::iterator edgeIt = myEdges.begin(); edgeIt
				!= myEdges.end(); ++edgeIt) {

			const AutomatonNode& succNodeObj_ref = edgeIt->getSuccessorNode();
			Node_ptr succNode = access_method::access(succNodeObj_ref);
			label_id_to_index_map_It = label_id_to_index_map.find(
					edgeIt->getEdgeLabelId());
			if (label_id_to_index_map_It != label_id_to_index_map.end())
				DFAmapIterator2->second.at(label_id_to_index_map_It->second)
						= succNode;
		}
	}

	Node_ptr secondDFAinitialNode = secondDFA_copy.getInitialNodes().at(0);
	std::set < Node_ptr > secondDFAacceptNodesSet;
	std::vector < Node_ptr > secondDFAacceptNodeList
			= secondDFA_copy.getAcceptNodes();
	for (std::vector<Node_ptr>::iterator it = secondDFAacceptNodeList.begin(); it
			!= secondDFAacceptNodeList.end(); ++it)
		secondDFAacceptNodesSet.insert(*it);

	/** Create new UnionDFA */
	DFAutomaton UnionDFA;
	std::queue < std::pair<Node_ptr, Node_ptr> > reachable_nodepairQueue;
	std::map < std::pair<Node_ptr, Node_ptr>, Node_ptr > NodeSet_to_Node_map;
	std::map<std::pair<Node_ptr, Node_ptr>, Node_ptr>::iterator
			NodeSet_to_Node_map_It1, NodeSet_to_Node_map_It2;

	std::pair < Node_ptr, Node_ptr > reachable_nodepair;
	reachable_nodepair = std::pair<Node_ptr, Node_ptr>(thisinitialNode,
			secondDFAinitialNode);

	const AutomatonNode& newNodeObj_ref = UnionDFA.createNode();
	Node_ptr newNode = access_method::access(newNodeObj_ref);

	UnionDFA.setInitialNode(*newNode);
	NodeSet_to_Node_map.insert(std::pair<std::pair<Node_ptr, Node_ptr>,
			Node_ptr>(reachable_nodepair, newNode));
	reachable_nodepairQueue.push(reachable_nodepair);

	/** flags to check if either one or both of the nodes in a pair are accept nodes */
	bool acceptNode1, acceptNode2;

	acceptNode1 = false;
	acceptNode2 = false;

	while (!reachable_nodepairQueue.empty()) {
		reachable_nodepair = reachable_nodepairQueue.front();
		Node_ptr firstNode, secondNode;
		firstNode = reachable_nodepair.first;
		secondNode = reachable_nodepair.second;
		reachable_nodepairQueue.pop();

		//!< node1 is an accepting node
		if (thisacceptNodesSet.count(firstNode) > 0)
			acceptNode1 = true;

		//!< node2 is an accepting node
		if (secondDFAacceptNodesSet.count(secondNode) > 0)
			acceptNode2 = true;

		/** Get the iterator */
		NodeSet_to_Node_map_It1 = NodeSet_to_Node_map.find(reachable_nodepair);

		if (NodeSet_to_Node_map_It1 == NodeSet_to_Node_map.end())
			throw std::runtime_error("@Union: NodeSet is not there in the map.");

		if (acceptNode1 || acceptNode2)
			UnionDFA.setAcceptNode(*(NodeSet_to_Node_map_It1->second));

		acceptNode1 = acceptNode2 = false;

		/** Get the vectors for both nodes in nodepair */
		DFAmapIterator1 = thisDFAmap.find(firstNode);
		DFAmapIterator2 = secondDFAmap.find(secondNode);

		std::set<label_cache::label_id_type>::size_type k = 0;
		for (std::set<label_cache::label_id_type>::iterator lIt =
				labelIds.begin(); lIt != labelIds.end(); ++lIt, ++k) {
			Node_ptr node1 = DFAmapIterator1->second.at(k);
			Node_ptr node2 = DFAmapIterator2->second.at(k);

			/** Create a new pair */
			reachable_nodepair = std::pair<Node_ptr, Node_ptr>(node1, node2);
			NodeSet_to_Node_map_It2 = NodeSet_to_Node_map.find(
					reachable_nodepair);

			/** If it is there in the map, add as a successor node */
			if (NodeSet_to_Node_map_It2 != NodeSet_to_Node_map.end()) {
				UnionDFA.addEdge(*(NodeSet_to_Node_map_It1->second),
						*(NodeSet_to_Node_map_It2->second), 1,
						label_cache::getLabel(*lIt));
			} //!< end-if

			/** else, create a new node corresponding to this pair,
			 * add to the map and finally add as a successor node */
			else {
				const AutomatonNode& newNodeObj_ref = UnionDFA.createNode();
				Node_ptr newNode = access_method::access(newNodeObj_ref);

				NodeSet_to_Node_map.insert(std::pair<std::pair<Node_ptr,
						Node_ptr>, Node_ptr>(reachable_nodepair, newNode));

				UnionDFA.addEdge(*(NodeSet_to_Node_map_It1->second), *newNode,
						1, label_cache::getLabel(*lIt));
				reachable_nodepairQueue.push(reachable_nodepair);
			} //!< end-else
		} //!< end-for (labelIds)
	}

	/** Set the alphabet */
	UnionDFA.setAlphabetIds(labelIds);

	/** Minimize Union DFA */
	DFAutomaton minUnion = UnionDFA.minimize();

	return minUnion;
}

DFAutomaton DFAutomaton::computeIntersection(const DFAutomaton& secondDFA) const {

	DFAutomaton firstDFA_copy(*this);
	DFAutomaton secondDFA_copy(secondDFA);

	std::set < label_cache::label_id_type > second_labelIdSet
			= secondDFA_copy.getAlphabetIds();
	std::set < label_cache::label_id_type > first_labelIdSet
			= firstDFA_copy.getAlphabetIds();

	firstDFA_copy.extend(second_labelIdSet);
	secondDFA_copy.extend(first_labelIdSet);

	first_labelIdSet.clear();
	second_labelIdSet.clear();

	/** Return the empty DFA as intersection */
	if (firstDFA_copy.getNodeList().empty())
		return firstDFA_copy;
	else if (secondDFA_copy.getNodeList().empty())
		return secondDFA_copy;

	/** label id's in an alphabet are mapped to their corresponding
	 * indices in the alphabet.	 */
	std::set < label_cache::label_id_type > labelIds
			= firstDFA_copy.getAlphabetIds();
	std::map<label_cache::label_id_type, unsigned int> label_id_to_index_map;
	std::map<label_cache::label_id_type, unsigned int>::iterator
			label_id_to_index_map_It;
	std::set<label_cache::label_id_type>::size_type i = 0;
	for (std::set<label_cache::label_id_type>::iterator lIt = labelIds.begin(); lIt
			!= labelIds.end(); ++lIt, ++i) {
		label_id_to_index_map.insert(std::pair<label_cache::label_id_type,
				unsigned int>(*lIt, i));
	}

	/** DFAmap - maps each DFA node to a vector of nodes.
	 *
	 *  For a node, entry at Index i in its corresponding vector
	 *  represents the successor node reached from transition on label
	 *  at index i in DFAalphabet. */
	std::map<Node_ptr, std::vector<Node_ptr> > thisDFAmap, secondDFAmap;
	std::map<Node_ptr, std::vector<Node_ptr> >::iterator DFAmapIterator1,
			DFAmapIterator2;

	std::vector < Node_ptr > thisNodeList = firstDFA_copy.getNodeList();

	for (std::vector<Node_ptr>::iterator it = thisNodeList.begin(); it
			!= thisNodeList.end(); ++it) {

		/** Initialize DFAmap for each node with an emptyDFANodeVector, which
		 *  will be filled later. */
		std::vector < Node_ptr > emptyDFANodeVector(labelIds.size());
		thisDFAmap.insert(std::pair<Node_ptr, std::vector<Node_ptr> >(*it,
				emptyDFANodeVector));

		DFAmapIterator1 = thisDFAmap.find(*it);

		if (DFAmapIterator1 == thisDFAmap.end())
			throw std::runtime_error(
					"@Intersect: Node is not there in thisDFAmap.");

		std::vector < AutomatonEdge > myEdges = (*it)->getEdges();

		/** fill emptyDFANodeVector */
		for (std::vector<AutomatonEdge>::iterator edgeIt = myEdges.begin(); edgeIt
				!= myEdges.end(); ++edgeIt) {

			const AutomatonNode& succNodeObj_ref = edgeIt->getSuccessorNode();
			Node_ptr succNode = access_method::access(succNodeObj_ref);
			label_id_to_index_map_It = label_id_to_index_map.find(
					edgeIt->getEdgeLabelId());
			if (label_id_to_index_map_It != label_id_to_index_map.end())
				DFAmapIterator1->second.at(label_id_to_index_map_It->second)
						= succNode;
		}
	}

	/** get initial node and accept nodes */
	Node_ptr thisinitialNode = firstDFA_copy.getInitialNodes().at(0);
	std::set < Node_ptr > thisacceptNodesSet;
	std::vector < Node_ptr > thisacceptNodeList
			= firstDFA_copy.getAcceptNodes();
	for (std::vector<Node_ptr>::iterator it = thisacceptNodeList.begin(); it
			!= thisacceptNodeList.end(); ++it)
		thisacceptNodesSet.insert(*it);

	/** Repeat the same for secondDFA */
	std::vector < Node_ptr > secondDFANodeList = secondDFA_copy.getNodeList();

	for (std::vector<Node_ptr>::iterator it = secondDFANodeList.begin(); it
			!= secondDFANodeList.end(); ++it) {

		/** Initialize DFAmap for each node with an emptyDFANodeVector to
		 *  be filled later. */
		std::vector < Node_ptr > emptyDFANodeVector(labelIds.size());
		secondDFAmap.insert(std::pair<Node_ptr, std::vector<Node_ptr> >(*it,
				emptyDFANodeVector));

		DFAmapIterator2 = secondDFAmap.find(*it);

		if (DFAmapIterator2 == secondDFAmap.end())
			throw std::runtime_error(
					"@Intersect: Node is not there in secondDFAmap.");

		std::vector < AutomatonEdge > myEdges = (*it)->getEdges();

		for (std::vector<AutomatonEdge>::iterator edgeIt = myEdges.begin(); edgeIt
				!= myEdges.end(); ++edgeIt) {

			const AutomatonNode& succNodeObj_ref = edgeIt->getSuccessorNode();
			Node_ptr succNode = access_method::access(succNodeObj_ref);
			label_id_to_index_map_It = label_id_to_index_map.find(
					edgeIt->getEdgeLabelId());
			if (label_id_to_index_map_It != label_id_to_index_map.end())
				DFAmapIterator2->second.at(label_id_to_index_map_It->second)
						= succNode;
		}
	}

	Node_ptr secondDFAinitialNode = secondDFA_copy.getInitialNodes().at(0);
	std::set < Node_ptr > secondDFAacceptNodesSet;
	std::vector < Node_ptr > secondDFAacceptNodeList
			= secondDFA_copy.getAcceptNodes();
	for (std::vector<Node_ptr>::iterator it = secondDFAacceptNodeList.begin(); it
			!= secondDFAacceptNodeList.end(); ++it)
		secondDFAacceptNodesSet.insert(*it);

	/** Create new IntersectDFA */
	DFAutomaton IntersectDFA;
	std::queue < std::pair<Node_ptr, Node_ptr> > reachable_nodepairQueue;
	std::map < std::pair<Node_ptr, Node_ptr>, Node_ptr > NodeSet_to_Node_map;
	std::map<std::pair<Node_ptr, Node_ptr>, Node_ptr>::iterator
			NodeSet_to_Node_map_It1, NodeSet_to_Node_map_It2;

	/** Set the alphabet */
	IntersectDFA.setAlphabetIds(labelIds);

	std::pair < Node_ptr, Node_ptr > reachable_nodepair;
	reachable_nodepair = std::pair<Node_ptr, Node_ptr>(thisinitialNode,
			secondDFAinitialNode);

	const AutomatonNode& newNodeObj_ref = IntersectDFA.createNode();
	Node_ptr newNode = access_method::access(newNodeObj_ref);

	IntersectDFA.setInitialNode(*newNode);
	NodeSet_to_Node_map.insert(std::pair<std::pair<Node_ptr, Node_ptr>,
			Node_ptr>(reachable_nodepair, newNode));
	reachable_nodepairQueue.push(reachable_nodepair);

	/** flags to check if either one or both of the nodes in a pair are accept nodes */
	bool acceptNode1, acceptNode2;

	acceptNode1 = false;
	acceptNode2 = false;

	while (!reachable_nodepairQueue.empty()) {
		reachable_nodepair = reachable_nodepairQueue.front();
		Node_ptr firstNode, secondNode;
		firstNode = reachable_nodepair.first;
		secondNode = reachable_nodepair.second;
		reachable_nodepairQueue.pop();

		//!< node1 is an accepting node
		if (thisacceptNodesSet.count(firstNode) > 0)
			acceptNode1 = true;

		//!< node2 is an accepting node
		if (secondDFAacceptNodesSet.count(secondNode) > 0)
			acceptNode2 = true;

		/** Get the iterator */
		NodeSet_to_Node_map_It1 = NodeSet_to_Node_map.find(reachable_nodepair);

		if (NodeSet_to_Node_map_It1 == NodeSet_to_Node_map.end())
			throw std::runtime_error(
					"@Intersect: NodeSet is not there in the map.");

		if (acceptNode1 && acceptNode2)
			IntersectDFA.setAcceptNode(*(NodeSet_to_Node_map_It1->second));

		acceptNode1 = acceptNode2 = false;

		/** Get the vectors for both nodes in nodepair */
		DFAmapIterator1 = thisDFAmap.find(firstNode);
		DFAmapIterator2 = secondDFAmap.find(secondNode);

		std::set<label_cache::label_id_type>::size_type k = 0;
		for (std::set<label_cache::label_id_type>::iterator lIt =
				labelIds.begin(); lIt != labelIds.end(); ++lIt, ++k) {
			Node_ptr node1 = DFAmapIterator1->second.at(k);
			Node_ptr node2 = DFAmapIterator2->second.at(k);

			/** Create a new pair */
			reachable_nodepair = std::pair<Node_ptr, Node_ptr>(node1, node2);
			NodeSet_to_Node_map_It2 = NodeSet_to_Node_map.find(
					reachable_nodepair);

			/** If it is there in the map, add as a successor node */
			if (NodeSet_to_Node_map_It2 != NodeSet_to_Node_map.end()) {
				(NodeSet_to_Node_map_It1->second)->addSuccNode(
						*(NodeSet_to_Node_map_It2->second), 1, *lIt);
			} //!< end-if

			/** else, create a new node corresponding to this pair
			 * add to the map and finally add as a successor node */
			else {
				const AutomatonNode& newNodeObj_ref = IntersectDFA.createNode();
				Node_ptr newNode = access_method::access(newNodeObj_ref);

				NodeSet_to_Node_map.insert(std::pair<std::pair<Node_ptr,
						Node_ptr>, Node_ptr>(reachable_nodepair, newNode));

				(NodeSet_to_Node_map_It1->second)->addSuccNode(*newNode, 1,
						*lIt);
				reachable_nodepairQueue.push(reachable_nodepair);
			} //!< end-else
		} //!< end-for (labelIds)
	}

	DFAutomaton minIntersect;
	/** Minimize Intersection DFA */
	if (!(IntersectDFA.is_empty()))
		minIntersect = IntersectDFA.minimize();

	return minIntersect;
}

DFAutomaton DFAutomaton::minimize() const {

	DFAutomaton thisDFA_copy(*this);

	std::vector < Node_ptr > NodeList = thisDFA_copy.getNodeList();

	if (NodeList.empty())
		throw std::runtime_error("@minimize - DFA is empty.");
	else
		NodeList.clear();

	/** Make this DFA complete, if not already. */
	thisDFA_copy.makeComplete();

	NodeList = thisDFA_copy.getNodeList();

	/** DFAmap - maps each DFA node to a vector of nodes.
	 *
	 *  For a node, entry at Index i in its corresponding vector
	 *  represents the successor node reached from transition on label
	 *  at index i in DFAalphabet.
	 */
	std::map < Node_ptr, std::vector<Node_ptr> > DFAmap;
	std::map<Node_ptr, std::vector<Node_ptr> >::iterator DFAmapIterator1,
			DFAmapIterator2;

	std::vector < Node_ptr > accept_nodes_list = thisDFA_copy.getAcceptNodes();

	if (accept_nodes_list.size() == 0)
		throw std::runtime_error("@minimize - DFA has no accept state.");

	/** label id's in an alphabet are mapped to their corresponding
	 * indices in the alphabet.
	 *
	 * @note Useful because we are dealing with a complete DFA.
	 * each node has a vector where an entry at index i depicts
	 * its successor node reachable on label at index i in its alphabet.
	 */
	std::set < label_cache::label_id_type > labelIds
			= thisDFA_copy.getAlphabetIds();
	std::map<label_cache::label_id_type, unsigned int> label_id_to_index_map;
	std::map<label_cache::label_id_type, unsigned int>::iterator
			label_id_to_index_mapIt;
	std::set<label_cache::label_id_type>::size_type i = 0;
	for (std::set<label_cache::label_id_type>::iterator lIt = labelIds.begin(); lIt
			!= labelIds.end(); ++lIt, ++i)
		label_id_to_index_map.insert(std::pair<label_cache::label_id_type,
				unsigned int>(*lIt, i));

	Node_ptr initNode;
	std::vector < Node_ptr > reachable_nodes;

	if (this->getInitialNodes().size() == 1) {
		initNode = thisDFA_copy.getInitialNodes().at(0);

		initNode->setStatus(AutomatonNode::IS_VISITED);
		reachable_nodes.push_back(initNode); // reachable_nodes is initialized with initNode.

	} else
		throw std::runtime_error("@minimize - DFA either has no initial "
			"state, or more than one initial state");

	/** Representative nodes for respective sets */
	Node_ptr accept_node_set_rep, non_accept_node_set_rep;

	/** Keep track of accepting and non-accepting nodes.
	 *
	 * @note These structures are updated while computing reachable
	 *  set of nodes, and used later for computation of equivalent
	 *  states' sets. */
	std::set<Node_ptr> accept_node_set, non_accept_node_set;

	/** maps a state to its representative state */
	std::map < Node_ptr, Node_ptr > rep_state_map;
	std::map<Node_ptr, Node_ptr>::iterator rep_state_mapIterator;

	/** maps a node to its corresponding equivalent state set id.
	 *
	 * @note Used while computing Equivalent states' sets.
	 * Initially, all accept nodes are assigned id 1.
	 * non-accept_nodes are assigned id 0. */
	std::map<Node_ptr, unsigned int> node_to_set_id_map;

	/** list to set conversion
	 * @note To avoid scan of whole accept_nodes_list every time
	 * a node is checked for an accept node. */
	for (std::vector<Node_ptr>::iterator aNodeIt = accept_nodes_list.begin(); aNodeIt
			!= accept_nodes_list.end(); ++aNodeIt)
		accept_node_set.insert(*aNodeIt);

	accept_nodes_list.clear();

	accept_node_set_rep = *(accept_node_set.begin());

	/** for each reachable node */
	for (std::vector<Node_ptr>::size_type rsize = 0; rsize
			< reachable_nodes.size(); ++rsize) {
		Node_ptr rNode = reachable_nodes.at(rsize);

		/** Each node in rep_state_map is mapped to first element of its
		 * corresponding set. */
		if (accept_node_set.count(rNode) > 0) {
			rep_state_map.insert(std::pair<Node_ptr, Node_ptr>(rNode,
					accept_node_set_rep));
			node_to_set_id_map.insert(std::pair<Node_ptr, unsigned int>(rNode,
					1));

		} else {
			non_accept_node_set.insert(rNode);
			if (non_accept_node_set.size() == 1)
				non_accept_node_set_rep = *(non_accept_node_set.begin());
			rep_state_map.insert(std::pair<Node_ptr, Node_ptr>(rNode,
					non_accept_node_set_rep));
			node_to_set_id_map.insert(std::pair<Node_ptr, unsigned int>(rNode,
					0));
		}

		std::vector < AutomatonEdge > myEdges = rNode->getEdges();

		/** Initialize DFAmap for each node with an emptyDFANodeVector, which
		 * will be filled later. */
		std::vector < Node_ptr > emptyDFANodeVector(labelIds.size());
		DFAmap.insert(std::pair<Node_ptr, std::vector<Node_ptr> >(rNode,
				emptyDFANodeVector));
		DFAmapIterator1 = DFAmap.find(rNode);

		for (std::vector<AutomatonEdge>::iterator edgeIt = myEdges.begin(); edgeIt
				!= myEdges.end(); ++edgeIt) {

			/** fill DFANodeVector */
			const AutomatonNode& succNodeObj_ref = edgeIt->getSuccessorNode();
			Node_ptr succNode = access_method::access(succNodeObj_ref);
			label_id_to_index_mapIt = label_id_to_index_map.find(
					edgeIt->getEdgeLabelId());
			if (label_id_to_index_mapIt != label_id_to_index_map.end())
				DFAmapIterator1->second.at(label_id_to_index_mapIt->second)
						= succNode;
			else
				throw std::runtime_error(
						"@minimize: label is not there in the map.");

			/** mark adjacentNode as reachable. */
			if (succNode->getStatus() == AutomatonNode::IS_NOT_VISITED) {
				succNode->setStatus(AutomatonNode::IS_VISITED);
				reachable_nodes.push_back(succNode);
			}
		} //!< end-for (myEdges)

	} //!< end-for (reachable_nodes)

	std::set < std::set<Node_ptr> > equivalent_state_set;

	/** Only when both sets are non-empty, equivalent_state_set may get altered
	 * during the computation. */
	if (!non_accept_node_set.empty() && !accept_node_set.empty()) {

		/** Initially, equivalent_state_set consists of accept_nodes' set and
		 * non_accept_nodes' set. */
		equivalent_state_set.insert(non_accept_node_set);
		equivalent_state_set.insert(accept_node_set);

		/** Initially, map_id 0 is for non_accept_nodes and 1 for accept_nodes. */
		int map_id = 1;

		/** flag to check if equivalent_state_set modified during last iteration. */
		bool change = true;

		/** Until a fixed point is reached. */
		while (change) {

			change = false;

			std::set<std::set<Node_ptr> >::iterator eq_st_it =
					equivalent_state_set.begin();

			/** for each equivalent state set */
			for (; eq_st_it != equivalent_state_set.end();) {

				/** modified_eq_set is the modified eq_set from which
				 *  states are removed and added to new_sep_eq_state_set.
				 *
				 *  Initially old_eq_set and modified_eq_set are similar.
				 *  old_eq_set is only traversed to check if it needs
				 *  to be changed. Based on this check, modified_eq_set is
				 *  being changed during the implementation */
				std::set<Node_ptr> modified_eq_set, old_eq_set;
				modified_eq_set = old_eq_set = *eq_st_it;
				Node_ptr rep_old_eq_set; //!< representative state for old set

				/** Those sets with cardinality > 1 are considered
				 * since such sets can only be partitioned further. */
				if (old_eq_set.size() > 1) {

					/** set representative state for old set */
					rep_old_eq_set = *old_eq_set.begin();

					std::set<Node_ptr>::iterator state1, state2;
					state1 = state2 = old_eq_set.begin();

					rep_state_mapIterator = rep_state_map.find(*state1);
					if (rep_state_mapIterator != rep_state_map.end())
						rep_state_mapIterator->second = rep_old_eq_set;

					/** For each pair (state1,state2) - check if all their corresponding
					 * transitions lead to same equivalent state set.
					 *
					 * States from old_eq_set are removed and are added to new_sep_eq_set. */
					std::set < Node_ptr > new_sep_eq_set;
					Node_ptr rep_new_sep_eq_set; //!< representative state for newly separated set

					state2++; //!< iterator to the successor state of state1, to form a pair.
					map_id++; //!< each node in new_sep_eq_set will be assigned this new map_id.

					DFAmapIterator1 = DFAmap.find(*state1);

					/** NodeVector for state1 */
					std::vector < Node_ptr > myNodeVector1
							= DFAmapIterator1->second;

					for (; state2 != old_eq_set.end(); ++state2) {

						rep_state_mapIterator = rep_state_map.find(*state2);
						if (rep_state_mapIterator != rep_state_map.end())
							rep_state_mapIterator->second = rep_old_eq_set;

						DFAmapIterator2 = DFAmap.find(*state2);

						/** NodeVector for state2 */
						std::vector < Node_ptr > myNodeVector2
								= DFAmapIterator2->second;

						/** maps to check the set ids of adjacent nodes. */
						std::map<Node_ptr, unsigned int>::iterator temp_map_it,
								temp_map_it1, temp_map_it2;

						for (std::vector<Node_ptr>::size_type i = 0; i
								< myNodeVector1.size(); i++) {

							/** get the ids of the nodes at index i in each vector
							 * of state1 and state2 */
							temp_map_it1 = node_to_set_id_map.find(
									myNodeVector1.at(i));
							temp_map_it2 = node_to_set_id_map.find(
									myNodeVector2.at(i));

							/** if transitions on same label for state1 and state2
							 * lead to states with different ids i.e., two different sets,
							 * then split modified_eq_set and add separated out state to
							 * new_sep_eq_set. */

							if (temp_map_it1->second != temp_map_it2->second) {

								Node_ptr unaryNode = *state2;
								new_sep_eq_set.insert(unaryNode);
								modified_eq_set.erase(*state2);

								if (new_sep_eq_set.size() == 1)
									rep_new_sep_eq_set
											= *new_sep_eq_set.begin();

								/** set representative state for newly separated set */
								rep_state_mapIterator = rep_state_map.find(
										*state2);
								if (rep_state_mapIterator
										!= rep_state_map.end())
									rep_state_mapIterator->second
											= rep_new_sep_eq_set;
								break;
							} //!< if (temp_map_it1->second != temp_map_it2->second)
						} //!< end-for (myNodeVector1)
					} //!< end-for (state2)

					/** Only if a state is separated
					 *
					 * change is set true only if old_eq_set is modified.
					 * If so, equivalent_state_set has to be changed.
					 *
					 * First remove old_eq_set from equivalent_state_set
					 * and then add modified_eq_set to equivalent_state_set.
					 */
					if (new_sep_eq_set.size() > 0) {
						equivalent_state_set.insert(new_sep_eq_set);

						for (std::set<Node_ptr>::iterator sep_it =
								new_sep_eq_set.begin(); sep_it
								!= new_sep_eq_set.end(); ++sep_it) {
							std::map<Node_ptr, unsigned int>::iterator
									temp_map_it;
							temp_map_it = node_to_set_id_map.find(*sep_it);
							temp_map_it->second = map_id;
						}
						change = true;
						equivalent_state_set.erase(old_eq_set);
						equivalent_state_set.insert(modified_eq_set);

						eq_st_it = equivalent_state_set.begin(); // reset the iterator for next iteration
					} else
						eq_st_it++;
				} //!<if (old_eq_set.size() > 1)

				/** else, every unary set has its own node as representative node */
				else {
					Node_ptr unaryNode = *(old_eq_set.begin());
					rep_state_mapIterator = rep_state_map.find(unaryNode);
					if (rep_state_mapIterator->second != unaryNode)
						rep_state_mapIterator->second = unaryNode;
					eq_st_it++;
				}
			} //!< end-for (equivalent_state_set)
		} //!< end-while
	} else {
		//!< std::cout << std::endl << "@minimize Either One or both sets are empty";
	}

	/** Create a new DFA */
	DFAutomaton minDFA;

	/** final node list and nodeIds */
	std::vector < Node_ptr > final_list;
	std::vector < new_node_id::node_id_type > final_nodeIds;

	for (std::vector<Node_ptr>::size_type i = 0; i < reachable_nodes.size(); ++i) {

		Node_ptr rNode = reachable_nodes.at(i);
		rNode->setStatus(AutomatonNode::IS_NOT_VISITED); //!< Reset the flag.
		rep_state_mapIterator = rep_state_map.find(rNode);
		if (rep_state_mapIterator != rep_state_map.end()) {

			bool add_this_node = false;

			/** Only representative nodes are added as each equivalent state set
			 *  has just one rep node, and rest of the nodes in that set
			 *  are replaced with this particular node. */

			if (rNode == rep_state_mapIterator->second)
				add_this_node = true;

			/** Check for initial node */
			if (rNode == initNode) {
				minDFA.setInitialNode(*(rep_state_mapIterator->second));
			}

			/** Check for final nodes */
			if (accept_node_set.count(rNode) > 0) {
				if (rNode == rep_state_mapIterator->second) {
					minDFA.setAcceptNode(*rNode);
				}
			}
			/** Replace each adjacent Node by its representative node, if both are different. */
			std::vector < AutomatonEdge > myEdges = rNode->getEdges();
			for (std::vector<AutomatonEdge>::size_type j = 0; j
					< myEdges.size(); ++j) {
				const AutomatonNode& succNodeObj_ref =
						myEdges[j].getSuccessorNode();
				Node_ptr succNode = access_method::access(succNodeObj_ref);
				rep_state_mapIterator = rep_state_map.find(succNode);
				if (rep_state_mapIterator->second != succNode)
					(rNode)->setSuccNodeAtIdx(*(rep_state_mapIterator->second),
							j);
				else {
					//!< No need to update adjacent Node as it is the representative node for its set.
				}
			}
			if (add_this_node)
				final_list.push_back(rNode);
		}
	}

	std::vector < new_node_id::node_id_type > NodeIds
			= thisDFA_copy.getNodeIds();
	NodeIds.resize(final_list.size());

	/** Only if it has at least one accepting state */
	if (minDFA.getAcceptNodes().size() > 0) {

		/** Update node list. */
		minDFA.resetNodeList(final_list);
		minDFA.resetNodeIds(NodeIds);
		for (std::set<label_cache::label_id_type>::iterator lIt =
				labelIds.begin(); lIt != labelIds.end(); ++lIt)
			minDFA.addToAlphabet(label_cache::getLabel(*lIt));
	} else
		throw std::runtime_error(
				"@minimize: minimized DFA doesn't have an accepting state.");

	/** minimize DFA */
	if (!(minDFA.is_empty()))
		return minDFA;
	else
		throw std::runtime_error("@minimize: minimized DFA is empty.");
}

std::vector<std::string> DFAutomaton::getShortestWord() const {

	Node_ptr source = this->getInitialNodes().at(0);
	std::vector < Node_ptr > myAcceptNodes = this->getAcceptNodes();
	std::set < Node_ptr > myAcceptNodesSet;
	for (std::vector<Node_ptr>::iterator aIt = myAcceptNodes.begin(); aIt
			!= myAcceptNodes.end(); ++aIt)
		myAcceptNodesSet.insert(*aIt);

	std::map < Node_ptr, temp_info_type > temp_info_map;
	std::map<Node_ptr, temp_info_type>::iterator temp_info_mapIt;

	std::queue < Node_ptr > reachable_nodes_queue;
	source->setStatus(AutomatonNode::IS_VISITED);
	reachable_nodes_queue.push(source);
	temp_info_type temp_info;
	temp_info_map.insert(std::pair<Node_ptr, temp_info_type>(source, temp_info));
	bool found = false;
	Node_ptr shortest_accepting_node;

	while (!reachable_nodes_queue.empty()) {

		Node_ptr rNode = reachable_nodes_queue.front();
		reachable_nodes_queue.pop();

		/** If an accepting node is found, break the loop. */
		if (myAcceptNodesSet.count(rNode) > 0) {
			found = true;
			shortest_accepting_node = rNode;
			break;
		}

		std::vector < AutomatonEdge > myEdges = rNode->getEdges();

		for (std::vector<AutomatonEdge>::iterator eIt = myEdges.begin(); eIt
				!= myEdges.end(); ++eIt) {

			distance_type uv_weight = (*eIt).getEdgeCost();
			const AutomatonNode& succNodeObj_ref = (*eIt).getSuccessorNode();
			Node_ptr mySuccNode = access_method::access(succNodeObj_ref);
			label_cache::label_id_type lId = (*eIt).getEdgeLabelId();

			if (mySuccNode->getStatus() == AutomatonNode::IS_NOT_VISITED) {
				mySuccNode->setStatus(AutomatonNode::IS_VISITED);
				temp_info_type temp_info;
				temp_info.parent = rNode;
				temp_info.label_id = lId;
				temp_info.distance = uv_weight; // Not relevant here.

				temp_info_map.insert(std::pair<Node_ptr, temp_info_type>(
						mySuccNode, temp_info));
				reachable_nodes_queue.push(mySuccNode);
			}
		}
	}

	/** Reset Nodes statuses */
	temp_info_mapIt = temp_info_map.begin();
	for (; temp_info_mapIt != temp_info_map.end(); ++temp_info_mapIt) {
		Node_ptr node = temp_info_mapIt->first;
		node->setStatus(AutomatonNode::IS_NOT_VISITED);
	}

	std::vector < std::string > shortest_word;
	if (found == false) {
		throw std::runtime_error(
				"@ShortestWord: There is no accepting state in the DFA.");
	} else {

		Node_ptr parent = shortest_accepting_node;
		std::vector < label_cache::label_id_type > Rev_shortest_word;
		temp_info_type temp_info;

		temp_info_mapIt = temp_info_map.find(parent);
		if (temp_info_mapIt != temp_info_map.end()) {
			temp_info = temp_info_mapIt->second;
			parent = temp_info.parent;
			label_cache::label_id_type label_id = temp_info.label_id;
			Rev_shortest_word.push_back(label_id);
		}

		while (parent != source) {
			temp_info_mapIt = temp_info_map.find(parent);
			if (temp_info_mapIt != temp_info_map.end()) {
				temp_info = temp_info_mapIt->second;
				parent = temp_info.parent;
				label_cache::label_id_type label_id = temp_info.label_id;
				Rev_shortest_word.push_back(label_id);
			}
		}
		std::vector<label_cache::label_id_type>::reverse_iterator rIt;
		for (rIt = Rev_shortest_word.rbegin(); rIt < Rev_shortest_word.rend(); ++rIt) {
			std::string label = label_cache::getLabel(*rIt);
			shortest_word.push_back(label);
		}
	}
	return shortest_word;
}

}

#endif /* DFAUTOMATON_HPP_ */
