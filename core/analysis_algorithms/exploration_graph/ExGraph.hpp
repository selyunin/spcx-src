/*
 * ExGraph.hpp
 *
 *  Created on: Apr 17, 2012
 *      Author: manish
 */

#ifndef EXGRAPH_HPP_
#define EXGRAPH_HPP_

#include "ExGraph.h"

namespace exploration_graph {

ExGraph::ExGraph() {
}

ExGraph::~ExGraph() {
}

ExGraph::ExGraph(const ExGraph& exGraph) {

	std::vector<ExNode_ptr> NodeList = exGraph.getNodeList();
	std::map<ExNode_ptr, ExNode_ptr> Nodes_map;
	std::map<ExNode_ptr, ExNode_ptr>::iterator Nodes_mapIt;
	for (std::vector<ExNode_ptr>::const_iterator nIt = NodeList.begin();
			nIt != NodeList.end(); ++nIt) {
		const ExGraphNode& newNodeObj_ref = createNode();
		ExNode_ptr exNode_shared_ptr = access_exNode_method::access(
				newNodeObj_ref);
		Nodes_map.insert(
				std::pair<ExNode_ptr, ExNode_ptr>(*nIt, exNode_shared_ptr));
	}

	for (std::vector<ExNode_ptr>::const_iterator nIt = NodeList.begin();
			nIt != NodeList.end(); ++nIt) {

		ExNode_ptr sNode;
		Nodes_mapIt = Nodes_map.find(*nIt);
		if (Nodes_mapIt != Nodes_map.end())
			sNode = Nodes_mapIt->second;
		std::vector<ExGraphTransition> myTransitions =
				(*nIt)->getTransitions();
		for (std::vector<ExGraphTransition>::iterator tIt =
				myTransitions.begin(); tIt != myTransitions.end(); ++tIt) {

			ExGraphTransition_Info_ptr myInfo = (*tIt).getInfo();
			ExNode_ptr mySuccNode = access_exNode_method::access(
					(*tIt).getSuccessorNode());
			trans_data_type tData = myInfo->getData();
			trans_type tType = myInfo->getType();

			ExNode_ptr tNode;
			Nodes_mapIt = Nodes_map.find(mySuccNode);
			if (Nodes_mapIt != Nodes_map.end())
				tNode = Nodes_mapIt->second;
			sNode->addTransition(*tNode, tType, tData);
		}
	}

	std::vector<ExNode_ptr> acceptNodes = exGraph.getAcceptNodes();
	for (std::vector<ExNode_ptr>::const_iterator nIt = acceptNodes.begin();
			nIt != acceptNodes.end(); ++nIt) {
		Nodes_mapIt = Nodes_map.find(*nIt);
		if (Nodes_mapIt != Nodes_map.end())
			setAcceptNode(*(Nodes_mapIt->second));
	}
	std::vector<ExNode_ptr> initialNodes = exGraph.getInitialNodes();
	for (std::vector<ExNode_ptr>::const_iterator nIt = initialNodes.begin();
			nIt != initialNodes.end(); ++nIt) {
		Nodes_mapIt = Nodes_map.find(*nIt);
		if (Nodes_mapIt != Nodes_map.end())
			setInitialNode(*(Nodes_mapIt->second));
	}
}

ExGraph ExGraph::computeReverse(bool only_visited_nodes) const {

	ExGraph reverseEG;

	std::map<ExNode_ptr, ExNode_ptr> Nodes_map;
	std::map<ExNode_ptr, ExNode_ptr>::iterator Nodes_mapIt;
	for (std::vector<ExNode_ptr>::const_iterator nIt = myNodeList.begin();
			nIt != myNodeList.end(); ++nIt) {

		if (only_visited_nodes == true) {
			if ((*nIt)->getStatus() == ExGraphNode::IS_VISITED) {
				const ExGraphNode& newNodeObj_ref = reverseEG.createNode();
				ExNode_ptr exNode_shared_ptr = access_exNode_method::access(
						newNodeObj_ref);
				Nodes_map.insert(
						std::pair<ExNode_ptr, ExNode_ptr>(*nIt,
								exNode_shared_ptr));
			}
		} else {
			const ExGraphNode& newNodeObj_ref = reverseEG.createNode();
			ExNode_ptr exNode_shared_ptr = access_exNode_method::access(
					newNodeObj_ref);
			Nodes_map.insert(
					std::pair<ExNode_ptr, ExNode_ptr>(*nIt, exNode_shared_ptr));

		}
	}

	for (std::vector<ExNode_ptr>::const_iterator nIt = myNodeList.begin();
			nIt != myNodeList.end(); ++nIt) {

		if (only_visited_nodes == true) {
			ExNode_ptr sNode;
			std::vector<ExGraphTransition> myTransitions;
			if ((*nIt)->getStatus() == ExGraphNode::IS_VISITED) {
				Nodes_mapIt = Nodes_map.find(*nIt);
				if (Nodes_mapIt != Nodes_map.end())
					sNode = Nodes_mapIt->second;
				myTransitions = (*nIt)->getTransitions();
			}
			for (std::vector<ExGraphTransition>::iterator tIt =
					myTransitions.begin(); tIt != myTransitions.end(); ++tIt) {

				ExGraphTransition_Info_ptr myInfo = (*tIt).getInfo();
				ExNode_ptr mySuccNode = access_exNode_method::access(
						(*tIt).getSuccessorNode());
				trans_data_type tData = myInfo->getData();
				trans_type tType = myInfo->getType();

				ExNode_ptr tNode;
				if (mySuccNode->getStatus() == ExGraphNode::IS_VISITED) {
					Nodes_mapIt = Nodes_map.find(mySuccNode);
					if (Nodes_mapIt != Nodes_map.end())
						tNode = Nodes_mapIt->second;
					tNode->addTransition(*sNode, tType, tData);
				}
			}
		} else {
			ExNode_ptr sNode;
			Nodes_mapIt = Nodes_map.find(*nIt);
			if (Nodes_mapIt != Nodes_map.end())
				sNode = Nodes_mapIt->second;
			std::vector<ExGraphTransition> myTransitions =
					(*nIt)->getTransitions();
			for (std::vector<ExGraphTransition>::iterator tIt =
					myTransitions.begin(); tIt != myTransitions.end(); ++tIt) {

				ExGraphTransition_Info_ptr myInfo = (*tIt).getInfo();
				ExNode_ptr mySuccNode = access_exNode_method::access(
						(*tIt).getSuccessorNode());
				trans_data_type tData = myInfo->getData();
				trans_type tType = myInfo->getType();

				ExNode_ptr tNode;
				Nodes_mapIt = Nodes_map.find(mySuccNode);
				if (Nodes_mapIt != Nodes_map.end())
					tNode = Nodes_mapIt->second;
				tNode->addTransition(*sNode, tType, tData);
			}
		}
	}

	for (std::vector<ExNode_ptr>::const_iterator nIt = myacceptNodes.begin();
			nIt != myacceptNodes.end(); ++nIt) {
		Nodes_mapIt = Nodes_map.find(*nIt);
		if (Nodes_mapIt != Nodes_map.end())
			reverseEG.setInitialNode(*(Nodes_mapIt->second));
	}

	for (std::vector<ExNode_ptr>::const_iterator nIt = myinitialNodes.begin();
			nIt != myinitialNodes.end(); ++nIt) {
		Nodes_mapIt = Nodes_map.find(*nIt);
		if (Nodes_mapIt != Nodes_map.end())
			reverseEG.setAcceptNode(*(Nodes_mapIt->second));
	}

	return reverseEG;
}

const ExGraphNode& ExGraph::createNode() {
	ExNode_ptr nNode(new ExGraphNode());
	node_id_type id = finite_automaton::new_node_id::create_id();
	myNodeList.push_back(nNode);
	myNodeIds.push_back(id);
	return *nNode;
}

void ExGraph::writeToDotFile(std::ofstream& outputFile) {

	std::map<ExNode_ptr, node_id_type> NodetoIdmap;
	std::map<ExNode_ptr, node_id_type>::iterator NodetoIdmapIt1, NodetoIdmapIt2;

	for (std::vector<ExNode_ptr>::size_type i = 0; i < myNodeList.size(); ++i)
		NodetoIdmap.insert(
				std::pair<ExNode_ptr, node_id_type>(myNodeList[i],
						myNodeIds[i]));

	std::set<ExNode_ptr> initialNodesSet;

	for (std::vector<ExNode_ptr>::iterator it = myinitialNodes.begin();
			it != myinitialNodes.end(); ++it)
		initialNodesSet.insert(*it);

	outputFile << "digraph language {" << std::endl;
	outputFile << "nodesep=.5;" << std::endl;
	outputFile << "rankdir = TB;" << std::endl;
	outputFile << "node [shape = doublecircle];";

	for (std::vector<ExNode_ptr>::iterator it = myacceptNodes.begin();
			it != myacceptNodes.end(); ++it) {
		NodetoIdmapIt1 = NodetoIdmap.find(*it);
		outputFile << NodetoIdmapIt1->second << " ";
	}
	outputFile << ";" << std::endl;

	outputFile << "node [shape = circle];" << std::endl;

	for (std::vector<ExNode_ptr>::iterator it = myinitialNodes.begin();
			it != myinitialNodes.end(); ++it) {
		NodetoIdmapIt1 = NodetoIdmap.find(*it);
		outputFile << "init";
		outputFile << " -> " << NodetoIdmapIt1->second << ";" << std::endl;
	}

	for (std::vector<ExNode_ptr>::size_type i = 0; i < myNodeList.size(); ++i) {

		NodetoIdmapIt1 = NodetoIdmap.find(myNodeList[i]);
		std::vector<ExGraphTransition> Transitions =
				myNodeList[i]->getTransitions();

		for (std::vector<ExGraphTransition>::size_type j = 0;
				j < Transitions.size(); ++j) {
			NodetoIdmapIt2 = NodetoIdmap.find(
					access_exNode_method::access(
							Transitions[j].getSuccessorNode()));
			if (NodetoIdmapIt2 != NodetoIdmap.end()) {

				std::stringstream ss;
				ExGraphTransition_Info_ptr myInfo = Transitions[j].getInfo();
				trans_type tType = myInfo->getType();
				trans_data_type tData = myInfo->getData();
				trans_bad_type bad = myInfo->getBadType();
				if (tType == ExGraphTransition_Info::DISCRETE) {
					if (transId_type * tId = boost::get<transId_type>(&tData))
						ss << *tId;
					outputFile << NodetoIdmapIt1->second;
					outputFile << " -> " << NodetoIdmapIt2->second;
					outputFile << " [ label = \"" << ss.str() << "\"";
					if (bad == ExGraphTransition_Info::IS_BAD)
						outputFile << ", color = red ];";
					else
						outputFile << "];";
				} else if (tType == ExGraphTransition_Info::TIME_ELAPSE) {
					if (interval_type * interval = boost::get<interval_type>(
							&tData))
						ss << *interval;
					outputFile << NodetoIdmapIt1->second;
					outputFile << " -> " << NodetoIdmapIt2->second;
					outputFile << " [ label = \"" << ss.str() << "\"";
					if (bad == ExGraphTransition_Info::IS_BAD)
						outputFile << ", color = red ];";
					else
						outputFile << "];";
				} else if (tType == ExGraphTransition_Info::CONTAINMENT) {
					outputFile << NodetoIdmapIt1->second;
					outputFile << " -> " << NodetoIdmapIt2->second;
					outputFile << "[ style = dotted";
					if (bad == ExGraphTransition_Info::IS_BAD)
						outputFile << ", color = red ];";
					else
						outputFile << "];";
				}
				outputFile << std::endl;
			}
		}
	}
	outputFile << "}";
	return;
}

void ExGraph::print() {
	std::map<ExNode_ptr, node_id_type> NodetoIdmap;
	std::map<ExNode_ptr, node_id_type>::iterator NodetoIdmapIt;

	std::set<ExNode_ptr> acceptNodesSet, initialNodesSet;
	for (std::vector<ExNode_ptr>::iterator it = myacceptNodes.begin();
			it != myacceptNodes.end(); ++it)
		acceptNodesSet.insert(*it);

	for (std::vector<ExNode_ptr>::iterator it = myinitialNodes.begin();
			it != myinitialNodes.end(); ++it)
		initialNodesSet.insert(*it);

	for (std::vector<ExNode_ptr>::size_type i = 0; i < myNodeList.size(); ++i)
		NodetoIdmap.insert(
				std::pair<ExNode_ptr, node_id_type>(myNodeList[i],
						myNodeIds[i]));

	for (std::vector<ExNode_ptr>::size_type i = 0; i < myNodeList.size(); ++i) {
		NodetoIdmapIt = NodetoIdmap.find(myNodeList[i]);
		std::cout << std::endl << NodetoIdmapIt->second;
		if (acceptNodesSet.count(myNodeList[i]) > 0)
			std::cout << "-- Accept Node";
		if (initialNodesSet.count(myNodeList[i]) > 0)
			std::cout << "-- Initial Node";

		std::vector<ExGraphTransition> Transitions =
				myNodeList[i]->getTransitions();
		for (std::vector<ExGraphTransition>::size_type j = 0;
				j < Transitions.size(); ++j) {
			std::stringstream ss;
			ExGraphTransition_Info_ptr myInfo = Transitions[j].getInfo();
			trans_type tType = myInfo->getType();
			trans_data_type tData = myInfo->getData();
			trans_bad_type bad = myInfo->getBadType();
			if (tType == ExGraphTransition_Info::DISCRETE) {
				ss << "(DISCRETE ";
				if (transId_type * tId = boost::get<transId_type>(&tData))
					ss << *tId;
			} else if (tType == ExGraphTransition_Info::TIME_ELAPSE) {
				ss << "(TIME_ELAPSE ";
				if (interval_type * interval = boost::get<interval_type>(
						&tData))
					ss << *interval;
			} else if (tType == ExGraphTransition_Info::CONTAINMENT) {
				ss << "(CONTAINMENT ";
				if (transId_type * tId = boost::get<transId_type>(&tData))
					ss << *tId;
				else if (interval_type * interval = boost::get<interval_type>(
						&tData))
					ss << *interval;
			}
			ss << ")";
			std::cout << std::endl << "--" << ss.str() << "-->";
			if (bad == ExGraphTransition_Info::IS_BAD) {
				std::cout << " bad -->";
			} else
				std::cout << " not bad -->";

			NodetoIdmapIt = NodetoIdmap.find(
					access_exNode_method::access(
							Transitions[j].getSuccessorNode()));
			if (NodetoIdmapIt != NodetoIdmap.end()) {
				std::cout << NodetoIdmapIt->second;
			} else
				std::cout << "No adjacent Node";
		}
	}
	std::cout << std::endl;
	return;
}

void ExGraph::addContainmentTransition(const ExGraphNode& sNode, const ExGraphNode& tNode) {
	ExGraphNode& sNode_ref = const_cast<ExGraphNode&>(sNode);
	transId_type tId = transId_type(0);
	sNode_ref.addTransition(tNode, ExGraphTransition_Info::CONTAINMENT, tId);
}

void ExGraph::addDiscreteTransition(const ExGraphNode& sNode, const ExGraphNode& tNode,
		const transId_type& tId) {
	ExGraphNode& sNode_ref = const_cast<ExGraphNode&>(sNode);
	sNode_ref.addTransition(tNode, ExGraphTransition_Info::DISCRETE, tId);
}

void ExGraph::addTimeElapseTransition(const ExGraphNode& sNode, const ExGraphNode& tNode,
		const interval_type& tInterval) {
	ExGraphNode& sNode_ref = const_cast<ExGraphNode&>(sNode);
	sNode_ref.addTransition(tNode, ExGraphTransition_Info::TIME_ELAPSE, tInterval);
}

void ExGraph::setInitialNode(const ExGraphNode& iNode) {
	myinitialNodes.push_back(access_exNode_method::access(iNode));
}

void ExGraph::setAcceptNode(const ExGraphNode& aNode) {
	myacceptNodes.push_back(access_exNode_method::access(aNode));
}

void ExGraph::mergeInitialNodes() {

	const ExGraphNode& newNode = this->createNode();
	ExGraphNode& newNode_ref = const_cast<ExGraphNode&>(newNode);

	transId_type tId_0 = transId_type(0);

	trans_type tType = ExGraphTransition_Info::DISCRETE;

	for (std::vector<ExNode_ptr>::const_iterator initIt =
			myinitialNodes.begin(); initIt != myinitialNodes.end(); ++initIt)
		newNode_ref.addTransition(*(*initIt), tType, tId_0);

	std::vector<ExNode_ptr> newInitialNodes;

	newInitialNodes.push_back(access_exNode_method::access(newNode));
	this->resetInitialNodes(newInitialNodes);

	return;
}

void ExGraph::backToOriginalInitialNodes() {

	ExNode_ptr newInitialNode = myinitialNodes.at(0);

	std::vector<ExNode_ptr> initialNodes;
	std::vector<ExGraphTransition> myTransitions =
			newInitialNode->getTransitions();

	for (std::vector<ExGraphTransition>::iterator tIt = myTransitions.begin();
			tIt != myTransitions.end(); ++tIt) {
		ExNode_ptr succNode = access_exNode_method::access(
				(*tIt).getSuccessorNode());
		initialNodes.push_back(succNode);
	}
	this->resetInitialNodes(initialNodes);
	myNodeList.erase(myNodeList.begin() + (myNodeList.size() - 1));
	myNodeIds.erase(myNodeIds.begin() + (myNodeIds.size() - 1));
	return;
}

ExGraph::shortest_path_info_type ExGraph::getShortestPath() {

	shortest_path_info_type shortest_path_info;

	this->mergeInitialNodes();

	ExNode_ptr myInitialNode = myinitialNodes.at(0);

	std::map<ExNode_ptr, temp_info_type> temp_info_map;
	std::map<ExNode_ptr, temp_info_type>::iterator temp_info_mapIt;

	distance_type max_distance = distance_type(1000);

	for (std::vector<ExNode_ptr>::iterator nIt = myNodeList.begin();
			nIt != myNodeList.end(); ++nIt) {
		temp_info_type temp_info;
		ExGraphTransition emptyTransition;

		if (*nIt == myInitialNode)
			temp_info.distance = distance_type(0);
		else
			temp_info.distance = max_distance;

		temp_info.transition = emptyTransition;

		temp_info_map.insert(
				std::pair<ExNode_ptr, temp_info_type>(*nIt, temp_info));
	}

	std::queue<ExNode_ptr> reachable_nodes;
	reachable_nodes.push(myInitialNode);

	while (!reachable_nodes.empty()) {
		ExNode_ptr rNode = reachable_nodes.front();
		reachable_nodes.pop();

		distance_type u_distance = distance_type(0);
		temp_info_mapIt = temp_info_map.find(rNode);
		if (temp_info_mapIt != temp_info_map.end()) {
			temp_info_type temp_info = temp_info_mapIt->second;
			u_distance = temp_info.distance;
		}

		std::vector<ExGraphTransition> myTransitions =
				rNode->getTransitions();

		for (std::vector<ExGraphTransition>::iterator tIt =
				myTransitions.begin(); tIt != myTransitions.end(); ++tIt) {

			distance_type uv_weight;
			distance_type v_distance = distance_type(0);

			ExGraphTransition_Info_ptr myInfo = (*tIt).getInfo(); //!< get transition information
			ExNode_ptr mySuccNode = access_exNode_method::access(
					(*tIt).getSuccessorNode());

			// Get ExTransition type
			trans_type tType = myInfo->getType();

			if (tType == ExGraphTransition_Info::DISCRETE)
				uv_weight = distance_type(1);
			else if (tType == ExGraphTransition_Info::TIME_ELAPSE)
				uv_weight = distance_type(1);
			else
				uv_weight = distance_type(0);

			temp_info_mapIt = temp_info_map.find(mySuccNode);
			if (temp_info_mapIt != temp_info_map.end()) {
				temp_info_type temp_info = temp_info_mapIt->second;
				v_distance = temp_info.distance;
			}

			if (u_distance + uv_weight < v_distance) {
				(temp_info_mapIt->second).distance = (u_distance + uv_weight);
				(temp_info_mapIt->second).parent = rNode;
				(temp_info_mapIt->second).transition = *tIt;
				reachable_nodes.push(mySuccNode);
			}
		}
	}

	distance_type shortest_distance = max_distance;
	ExNode_ptr shortest_accept_node;
	for (std::vector<ExNode_ptr>::iterator aIt = myacceptNodes.begin();
			aIt != myacceptNodes.end(); ++aIt) {
		distance_type temp_distance = distance_type(0);
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

	ExNode_ptr parent = shortest_accept_node;
	ExNode_ptr source = shortest_accept_node;
	std::vector<ExGraphTransition> Rev_shortest_Path;
	temp_info_type temp_info;

	temp_info_mapIt = temp_info_map.find(parent);
	if (temp_info_mapIt != temp_info_map.end()) {
		temp_info = temp_info_mapIt->second;
		parent = temp_info.parent;
	}
	while (parent != myInitialNode) {
		ExGraphTransition_Info_ptr myInfo = temp_info.transition.getInfo();
		trans_type tType = myInfo->getType();
		if (tType != ExGraphTransition_Info::CONTAINMENT)
			Rev_shortest_Path.push_back(temp_info.transition);
		source = parent;
		temp_info_mapIt = temp_info_map.find(parent);
		if (temp_info_mapIt != temp_info_map.end()) {
			temp_info = temp_info_mapIt->second;
			parent = temp_info.parent;
		}
	}
	std::vector<ExGraphTransition> shortest_Path = ReverseTranstionVector(
			Rev_shortest_Path);

	shortest_path_info.initialExNode = source;
	shortest_path_info.Transitions = shortest_Path;

	this->backToOriginalInitialNodes();
	return shortest_path_info;

}
}

#endif /* EXGRAPH_HPP_ */
