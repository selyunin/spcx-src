/*
 * ExGraph.h
 *
 *  Created on: Jan 5, 2012
 *      Author: goyal
 */

#ifndef EXGRAPH_H_
#define EXGRAPH_H_

#include "boost/shared_ptr.hpp"
#include <boost/enable_shared_from_this.hpp>
#include <sstream>
#include <stack>
#include <map>
#include <queue>
#include <stdexcept>
#include "ExGraphTransition.h"
#include "extern/finite_automaton/new_node_id.h"
#include "extern/finite_automaton/label_cache.h"
using namespace std;

namespace exploration_graph {

/** forward declarations */
typedef finite_automaton::new_node_id::node_id_type node_id_type;
typedef finite_automaton::label_cache::label_id_type label_id_type;

/** ExGraph class depicts an Exploration Graph which is consisted
 * of a set of nodes and a set of transitions. */
class ExGraph: public boost::enable_shared_from_this<ExGraph> {
public:

	/** Represents shortest path information in terms
	 * of an initial node and a list of transitions. */
	typedef struct {
		ExNode_ptr initialExNode;
		std::vector<ExGraphTransition> Transitions;
	} shortest_path_info_type;

	typedef unsigned int distance_type;

	typedef struct {
		ExNode_ptr parent;
		distance_type distance;
		ExGraphTransition transition;
	} temp_info_type;

	typedef boost::shared_ptr<ExGraph> ExGraph_ptr;

	/** Return a shared_ptr to *this. */
	ExGraph_ptr get_ptr() {
		return boost::enable_shared_from_this<ExGraph>::shared_from_this();
	}
	;

	/** default constructor */
	ExGraph();

	/** Construct a deep copy of exGraph.
	 *
	 * O(n+m)log n. */
	ExGraph(const ExGraph& exGraph);

	/** ExGraph destructor */
	virtual ~ExGraph();

	/** Returns true if the language is empty, else false */
	bool is_empty();

	/** Create a node and return its reference. */
	const ExGraphNode& createNode();

	/** Add Containment transition.
	 *
	 * By default, the transition data is set to transition id 0. */
	void addContainmentTransition(const ExGraphNode& sNode, const ExGraphNode& tNode);

	/** Add Discrete Transition. */
	void addDiscreteTransition(const ExGraphNode& sNode, const ExGraphNode& tNode,
			const transId_type& tId);

	/** Add Time Elapse Transition. */
	void addTimeElapseTransition(const ExGraphNode& sNode, const ExGraphNode& tNode,
			const interval_type& tInterval);

	/** Add Initial node to the ExGraph
	 *
	 * @note Assumed that iNode is already there in EG. */
	void setInitialNode(const ExGraphNode& iNode);

	/** Add Accept node to the ExGraph
	 *
	 * @note Assumed that aNode is already there in EG. */
	void setAcceptNode(const ExGraphNode& aNode);

	/** clear the contents of an ExGraph. */
	void clear();

	/** Print an ExGraph.
	 *
	 * O(n+m)log n. */
	void print();

	/** Write exploration graph to a dot file.
	 *
	 * O(n+m)log n. */
	void writeToDotFile(std::ofstream& outputFile);

	/** Return shortest path in terms of an initial node
	 * and a list of transitions.
	 *
	 * @note Containment transitions are not included in the
	 * computed path as they are assumed to be of cost 0.
	 * while Discrete and Time_Elapse transitions are
	 * considered to be of cost 1.
	 *
	 * O(n+m)log n. */
	shortest_path_info_type getShortestPath();

private:
	std::vector<ExNode_ptr> myNodeList; /*!< List of nodes */
	std::vector<ExNode_ptr> myacceptNodes; /*!< List of accepting nodes */
	std::vector<node_id_type> myNodeIds; /*!< Nodes id's */
	std::vector<ExNode_ptr> myinitialNodes; /*!< List of initial nodes */

	/** Reverse the Transition vector computed in
	 * getShortestPath method. */
	std::vector<ExGraphTransition> ReverseTranstionVector(
			std::vector<ExGraphTransition> TransitionVector) {

		std::vector<ExGraphTransition> shortestpath;
		if (TransitionVector.empty())
			return shortestpath;

		std::vector<ExGraphTransition>::iterator tIt = TransitionVector.end();
		tIt--;
		for (; tIt != TransitionVector.begin(); tIt--)
			shortestpath.push_back(*tIt);

		shortestpath.push_back(*tIt);

		return shortestpath;
	}
	;

	friend class EG_NFA_interface;

protected:

	/** Reset initial nodes.
	 *
	 * Used during MergeInitialNodes and backToOriginalInitialNodes
	 * methods. */
	void resetInitialNodes(std::vector<ExNode_ptr> newInitialNodes) {
		myinitialNodes.clear();
		for (std::vector<ExNode_ptr>::iterator initIt = newInitialNodes.begin();
				initIt != newInitialNodes.end(); ++initIt)
			myinitialNodes.push_back(*initIt);
	}

	/**
	 * MergeInitialNodes and backToOriginalInitialNodes methods are
	 * used while computing a shortest path in an exploration graph.
	 *
	 * Before computing a shortest path, discrete transitions with id '0'
	 * are added from new initial node to each old initial node; and, after
	 * computing the path, new initial node is deleted; thus, yielding back
	 * our original exploration graph configuration.
	 */
	void mergeInitialNodes();

	void backToOriginalInitialNodes();

	/** Compute the reverse of an EG.
	 *
	 * only_visited_nodes is a flag. If true, only visited nodes and
	 * the transitions between them are accessed.
	 *
	 * O(n+m)log n */
	ExGraph computeReverse(bool only_visited_nodes) const;

	/** Return the vector of nodes in an ExGraph */
	std::vector<ExNode_ptr> getNodeList() const {
		return myNodeList;
	}

	/** Return a vector of initial nodes */
	std::vector<ExNode_ptr> getInitialNodes() const {
		return myinitialNodes;
	}

	/** Return a vector of accept nodes */
	std::vector<ExNode_ptr> getAcceptNodes() const {
		return myacceptNodes;
	}

	/** Get NodeIds of the graph. */
	std::vector<node_id_type> getNodeIds() const {
		return myNodeIds;
	}
	/** Reset node list */
	void resetNodeList(const std::vector<ExNode_ptr> nodeList) {
		myNodeList.clear();
		myNodeList = nodeList;
	}

	/** Reset nodes id's */
	void resetNodeIds(const std::vector<node_id_type> nodeIds) {
		myNodeIds.clear();
		myNodeIds = nodeIds;
	}

};
}
#include "ExGraph.hpp"
#endif /* EXGRAPH_H_ */
