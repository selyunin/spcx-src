#ifndef AUTOMATON_HPP_
#define AUTOMATON_HPP_

#include "Automaton.h"

namespace finite_automaton {

/** distance type is same as of cost type. */
typedef AutomatonEdge::cost_type distance_type;

/** temporary information type used during the
 *  shortest path computation. */
typedef struct {
	Node_ptr parent;
	distance_type distance;
	label_cache::label_id_type label_id;
} temp_info_type;

Automaton::Automaton() {
}

Automaton::~Automaton() {
}

bool Automaton::is_empty() {

	bool is_empty = false;
	if (myNodeList.empty())
		is_empty = true;
	else if (myinitialNodes.empty())
		is_empty = true;
	else if (myacceptNodes.empty())
		is_empty = true;
	else {
		is_empty = true;
		std::set < Node_ptr > myacceptNodesSet;
		for (std::vector<Node_ptr>::iterator aIt = myacceptNodes.begin(); aIt
				!= myacceptNodes.end(); ++aIt)
			myacceptNodesSet.insert(*aIt);

		for (std::vector<Node_ptr>::iterator iIt = myinitialNodes.begin(); iIt
				!= myinitialNodes.end() && is_empty; iIt++) {
			if (myacceptNodesSet.count(*iIt) > 0) {
				is_empty = false;
				break;
			} //!< end-if(myacceptNodesSet.count(*iIt) > 0)
			else {
				std::queue < Node_ptr > reachableNodes;
				reachableNodes.push(*iIt);

				while (!reachableNodes.empty() && is_empty) {
					Node_ptr front = reachableNodes.front();
					front->setStatus(AutomatonNode::IS_VISITED);
					reachableNodes.pop();
					std::vector < AutomatonEdge > Edges = front->getEdges();
					for (std::vector<AutomatonEdge>::iterator eIt =
							Edges.begin(); eIt != Edges.end(); ++eIt) {
						const AutomatonNode& succNodeObj_ref =
								(*eIt).getSuccessorNode();
						Node_ptr succNode = access_method::access(
								succNodeObj_ref);
						if (succNode->getStatus()
								== AutomatonNode::IS_NOT_VISITED) {
							if (myacceptNodesSet.count(succNode) > 0) {
								is_empty = false;
								break;
							} //!< end-if(myacceptNodesSet.count(adjNode) > 0)
							reachableNodes.push(succNode);
						} //!< end-if(adjNode->getMyStatus() == AutomatonNode::IS_NOT_VISITED)
					} //!< end-for
				} //!< end-while
			} //!< end-else if(myacceptNodesSet.count(*iIt) > 0)
		} //!< end-for
	} //!< end-else

	/** Reset visited status of every node */
	for (std::vector<Node_ptr>::iterator nIt = myNodeList.begin(); nIt
			!= myNodeList.end(); ++nIt)
		(*nIt)->setStatus(AutomatonNode::IS_NOT_VISITED);

	if (is_empty)
		return true;
	else
		return false;
}

Automaton::Automaton(const Automaton& automaton) {

	std::vector < Node_ptr > NodeList = automaton.getNodeList();
	std::set < label_cache::label_id_type > alphabetIds
			= automaton.getAlphabetIds();
	myalphabetIds = alphabetIds;
	std::map < Node_ptr, Node_ptr > Nodes_map;
	std::map<Node_ptr, Node_ptr>::iterator Nodes_mapIt;
	for (std::vector<Node_ptr>::const_iterator nIt = NodeList.begin(); nIt
			!= NodeList.end(); ++nIt) {
		const AutomatonNode& newNodeObj_ref = createNode();
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

void Automaton::print() {

	std::map < Node_ptr, new_node_id::node_id_type > NodetoIdmap;
	std::map<Node_ptr, new_node_id::node_id_type>::iterator NodetoIdmapIt;

	std::set<Node_ptr> acceptNodesSet, initialNodesSet;
	for (std::vector<Node_ptr>::iterator it = myacceptNodes.begin(); it
			!= myacceptNodes.end(); ++it)
		acceptNodesSet.insert(*it);

	for (std::vector<Node_ptr>::iterator it = myinitialNodes.begin(); it
			!= myinitialNodes.end(); ++it)
		initialNodesSet.insert(*it);

	for (std::vector<Node_ptr>::size_type i = 0; i < myNodeList.size(); ++i)
		NodetoIdmap.insert(std::pair<Node_ptr, new_node_id::node_id_type>(
				myNodeList[i], myNodeIds[i]));

	for (std::vector<Node_ptr>::size_type i = 0; i < myNodeList.size(); ++i) {
		NodetoIdmapIt = NodetoIdmap.find(myNodeList[i]);
		std::cout << std::endl << NodetoIdmapIt->second;
		if (acceptNodesSet.count(myNodeList[i]) > 0)
			std::cout << "-- Accept Node";
		if (initialNodesSet.count(myNodeList[i]) > 0)
			std::cout << "-- Initial Node";

		std::vector < AutomatonEdge > Edges = myNodeList[i]->getEdges();
		for (std::vector<AutomatonEdge>::size_type j = 0; j < Edges.size(); ++j) {
			std::cout << std::endl << "--" << label_cache::getLabel(
					Edges[j].getEdgeLabelId()) << "-->";
			const AutomatonNode& succNodeObj_ref = Edges[j].getSuccessorNode();
			Node_ptr succNode = access_method::access(succNodeObj_ref);
			NodetoIdmapIt = NodetoIdmap.find(succNode);
			if (NodetoIdmapIt != NodetoIdmap.end()) {
				std::cout << NodetoIdmapIt->second;
			} else
				std::cout << "No adjacent Node";
		}
	}
	std::cout << std::endl;
	return;
}

void Automaton::writeToDotFile(std::ofstream& outputFile) const {

	std::map < Node_ptr, new_node_id::node_id_type > NodetoIdmap;
	std::map<Node_ptr, new_node_id::node_id_type>::iterator NodetoIdmapIt1,
			NodetoIdmapIt2;

	for (std::vector<Node_ptr>::size_type i = 0; i < myNodeList.size(); ++i)
		NodetoIdmap.insert(std::pair<Node_ptr, new_node_id::node_id_type>(
				myNodeList[i], myNodeIds[i]));

	outputFile << "digraph language {" << std::endl;
	outputFile << "rankdir = LR;" << std::endl;
	outputFile << "d2tdocpreamble = \"\\usetikzlibrary{automata}\";"
			<< std::endl;
	outputFile
			<< "d2tfigpreamble = \"\\tikzstyle{every state} = [draw=blue!50, very thick, fill=blue!10]\";"
			<< std::endl;
	outputFile << "node [style = \"state\"];" << std::endl;
	outputFile << "edge [lblstyle=\"auto\"];" << std::endl;

	for (std::vector<Node_ptr>::const_iterator it = myinitialNodes.begin(); it
			!= myinitialNodes.end(); ++it) {
		NodetoIdmapIt1 = NodetoIdmap.find(*it);
		outputFile << NodetoIdmapIt1->second << " [style=\"state,initial\"];"
				<< std::endl;
	}

	for (std::vector<Node_ptr>::size_type i = 0; i < myNodeList.size(); ++i) {
		std::map < new_node_id::node_id_type, std::vector<std::string>
				> NodetoLablesmap;
		std::map<new_node_id::node_id_type, std::vector<std::string> >::iterator
				NodetoLablesmapIt;

		NodetoIdmapIt1 = NodetoIdmap.find(myNodeList[i]);
		std::vector < AutomatonEdge > Edges = myNodeList[i]->getEdges();

		/** fill NodetoLablesmap */
		for (std::vector<AutomatonEdge>::size_type j = 0; j < Edges.size(); ++j) {
			const AutomatonNode& succNodeObj_ref = Edges[j].getSuccessorNode();
			Node_ptr succNode = access_method::access(succNodeObj_ref);
			NodetoIdmapIt2 = NodetoIdmap.find(succNode);
			if (NodetoIdmapIt2 != NodetoIdmap.end()) {
				NodetoLablesmapIt
						= NodetoLablesmap.find(NodetoIdmapIt2->second);
				if (NodetoLablesmapIt != NodetoLablesmap.end())
					NodetoLablesmapIt->second.push_back(label_cache::getLabel(
							Edges[j].getEdgeLabelId()));
				else {
					std::vector < std::string > Labels;
					Labels.push_back(label_cache::getLabel(
							Edges[j].getEdgeLabelId()));
					NodetoLablesmap.insert(std::pair<new_node_id::node_id_type,
							std::vector<std::string> >(NodetoIdmapIt2->second,
							Labels));
				}
			}
		}
		for (NodetoLablesmapIt = NodetoLablesmap.begin(); NodetoLablesmapIt
				!= NodetoLablesmap.end(); ++NodetoLablesmapIt) {
			outputFile << NodetoIdmapIt1->second;
			outputFile << " -> " << NodetoLablesmapIt->first;
			std::vector < std::string > Labels = NodetoLablesmapIt->second;
			outputFile << " [ label = \"";
			std::vector<std::string>::iterator lIt = Labels.begin();
			outputFile << *lIt;
			++lIt;
			for (; lIt != Labels.end(); ++lIt)
				outputFile << ", " << *lIt;
			outputFile << "";
			outputFile << "\"];";
			outputFile << std::endl;
		}
	} //!< end-for nodes

	for (std::vector<Node_ptr>::const_iterator it = myacceptNodes.begin(); it
			!= myacceptNodes.end(); ++it) {
		NodetoIdmapIt1 = NodetoIdmap.find(*it);
		outputFile << NodetoIdmapIt1->second << " [style=\"state,accepting\"];"
				<< std::endl;
	}
	outputFile << "}";
	return;
}

void Automaton::clear() {
	myNodeList.clear();
	myacceptNodes.clear();
	myinitialNodes.clear();
	myalphabetIds.clear();
	myNodeIds.clear();
}

const AutomatonNode& Automaton::createNode() {
	Node_ptr nNode(new AutomatonNode());
	new_node_id::node_id_type id = new_node_id::create_id();
	myNodeList.push_back(nNode);
	myNodeIds.push_back(id);
	return *nNode;
}

void Automaton::addEdge(const AutomatonNode& sNode, const AutomatonNode& tNode,
		const cost_type& cost, const std::string& label) {

	label_cache::label_id_type lId = label_cache::getLabelId(label);

	AutomatonNode& sNode_ref = const_cast<AutomatonNode&> (sNode);

	if (lId != -1) {
		if ((myalphabetIds.count(lId)) > 0)
			sNode_ref.addSuccNode(tNode, cost, lId);
		else
			throw std::runtime_error(label
					+ " is not in the automaton alphabet.");
	} else
		throw std::runtime_error(label + " is not a valid label.");

}

void Automaton::setInitialNode(const AutomatonNode& iNode) {
	myinitialNodes.push_back(access_method::access(iNode));
}

void Automaton::setAcceptNode(const AutomatonNode& aNode) {
	myacceptNodes.push_back(access_method::access(aNode));
}

void Automaton::addToAlphabet(std::string label) {

	label_cache::label_id_type labelId;

	labelId = label_cache::getLabelId(label);

	if (labelId == -1)
		myalphabetIds.insert(label_cache::addLabel(label));
	else
		myalphabetIds.insert(labelId);
}

void Automaton::setAlphabet(const std::set<std::string> alphabet) {
	myalphabetIds.clear();
	for (std::set<std::string>::iterator it = alphabet.begin(); it
			!= alphabet.end(); ++it)
		this->addToAlphabet(*it);
}

std::set<std::string> Automaton::getAlphabet() const {
	std::set < std::string > alphabet;
	for (std::set<label_cache::label_id_type>::iterator lIt =
			myalphabetIds.begin(); lIt != myalphabetIds.end(); ++lIt)
		alphabet.insert(label_cache::getLabel(*lIt));
	return alphabet;
}

std::vector<std::string> Automaton::getShortestPath() const {

	Node_ptr source = myinitialNodes.at(0);
	distance_type max_distance = distance_type(MAX_DISTANCE);

	std::map < Node_ptr, temp_info_type > temp_info_map;
	std::map<Node_ptr, temp_info_type>::iterator temp_info_mapIt;

	for (std::vector<Node_ptr>::const_iterator nIt = myNodeList.begin(); nIt
			!= myNodeList.end(); ++nIt) {

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
	for (std::vector<Node_ptr>::const_iterator aIt = myacceptNodes.begin(); aIt
			!= myacceptNodes.end(); ++aIt) {
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

		std::vector<label_cache::label_id_type>::reverse_iterator rIt;
		for (rIt = Rev_shortest_path.rbegin(); rIt < Rev_shortest_path.rend(); ++rIt) {
			std::string label = label_cache::getLabel(*rIt);
			shortest_path.push_back(label);
		}

	}
	return shortest_path;
}
}
#endif /* AUTOMATON_HPP_ */
