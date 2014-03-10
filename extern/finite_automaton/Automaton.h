#ifndef AUTOMATON_H_
#define AUTOMATON_H_

/*********************************************************************
 * Automaton.h														 *
 * 																	 *
 * Copyright (C) 2012 by Verimag Research Lab, Grenoble, France.  	 *
 *																	 *
 * This program is free software; you can redistribute it and/or  	 *
 * modify it under the terms of the GNU Lesser General Public     	 *
 * License as published by the Free Software Foundation; either   	 *
 * version 2.1 of the License, or (at your option) any later version.*
 * 												  					 *
 * This program is distributed in the hope that it will be useful,	 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 	 *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU *
 * Lesser General Public License for more details.			  	 	 *
 *																	 *
 * You should have received a copy of the GNU Lesser General Public	 *
 * License along with this library; if not, write to the 			 *
 * Free	Software Foundation, Inc., 									 *
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.	 *
 *																	 *
 * @authors: Manish Goyal (manish.goyal@imag.fr)					 *
 * 			 Goran Frehse (goran.frehse@imag.fr)					 *
 *********************************************************************/

#define MAX_DISTANCE 10000

#include "boost/shared_ptr.hpp"
#include <boost/enable_shared_from_this.hpp>
#include <stack>
#include <queue>
#include <set>
#include <stdexcept>
#include "new_node_id.h"
#include "label_cache.h"
#include "AutomatonEdge.h"
#include "core/analysis_algorithms/exploration_graph/EG_NFA_interface.h"

namespace exploration_graph {
class EG_NFA_interface;
}

namespace finite_automaton {

/** forward declarations */
typedef AutomatonEdge::cost_type cost_type;

/** label id is initialized with 1.
 *
 * label id 0 is reserved for "EPS" transitions. */
label_cache::label_id_type label_cache::label_id = label_cache::label_id_type(1);

/** highestid initialized to 0. */
new_node_id::node_id_type new_node_id::highestid = new_node_id::node_id_type(0);

/** initiliazing label cache. */
std::map<label_cache::label_id_type, std::string> label_cache::label_id_to_str_map =
		label_cache::init_label_cache();

/**
 * Automaton class depicts the finite automaton which is consisted
 * of a set of nodes and a set of transitions.
 *
 * @note In our implementation, Automaton is stored as a vector of
 * node-pointers where, each node has its own list of outgoing
 * transitions. An automaton also contains a list of initial nodes,
 * accepting nodes and transition labels.
 */

class Automaton: public boost::enable_shared_from_this<Automaton> {
public:
	typedef boost::shared_ptr<Automaton> Automaton_ptr;

	/** Return a shared_ptr to *this. */
	Automaton_ptr get_ptr() {
		return boost::enable_shared_from_this<Automaton>::shared_from_this();
	}
	;

	/** Automaton constructors. */
	Automaton();

	/** Construct a deep copy of automaton.
	 * O(n+m)log n.*/
	Automaton(const Automaton& automaton);

	/** Automaton destructor. */
	virtual ~Automaton();

	/** Return true if the language is empty, else false.
	 * O(m). */
	bool is_empty();

	/** Create a new node of the automaton and return its reference. */
	const AutomatonNode& createNode();

	/** Add an edge to the automaton. */
	void addEdge(const AutomatonNode& sNode, const AutomatonNode& tNode,
			const cost_type& cost, const std::string& label);

	/** Add label to the automaton alphabet. */
	void addToAlphabet(std::string label);

	/** Return automaton alphabet. */
	std::set<std::string> getAlphabet() const;

	/** Set automaton alphabet. */
	void setAlphabet(const std::set<std::string> alphabet);

	/** Set the node as initial node.
	 *
	 * @note Assumed that iNode is already there in an automaton. */
	void setInitialNode(const AutomatonNode& iNode);

	/** Set the node as accepting node.
	 *
	 * @note Assumed that aNode is already there in an automaton. */
	void setAcceptNode(const AutomatonNode& aNode);

	/** clear the contents of a Automaton. */
	void clear();

	/** Print the automaton */
	void print();

	/** Write automaton to a dot file. */
	void writeToDotFile(std::ofstream& outputFile) const;

	/** Return shortest path, computed with respect to edge costs,
	 * as a vector of labels.
	 *
	 * Algorithm used is Bellman-ford.
	 * O(n+m)log n. */
	virtual std::vector<std::string> getShortestPath() const;

private:

	std::vector<Node_ptr> myNodeList; /*!< list of nodes */
	std::vector<Node_ptr> myacceptNodes; /*!< list of accept nodes */
	std::vector<new_node_id::node_id_type> myNodeIds; /*!< node id's */
	std::vector<Node_ptr> myinitialNodes; /*!< list of initial nodes */
	std::set<label_cache::label_id_type> myalphabetIds; /*!< alphabet is stored as a list of labels id's */

	friend class exploration_graph::EG_NFA_interface;

protected:

	/** Both methods are used while converting an NFA into a DFA,
	 * and DFA minimization. */

	/** Reset the node list. */
	void resetNodeList(const std::vector<Node_ptr> nodeList) {
		myNodeList.clear();
		myNodeList = nodeList;
	}

	/** Reset node id's of the automaton. */
	void resetNodeIds(const std::vector<new_node_id::node_id_type> NodeIds) {
		myNodeIds.clear();
		myNodeIds = NodeIds;
	}

	/** Reset initial node.
	 *
	 * @note Used while merging and separating multiple initial nodes
	 * in an NFA. */
	void resetInitialNode(const Node_ptr iNode) {
		myinitialNodes.clear();
		myinitialNodes.push_back(iNode);
	}

	/** Return a vector of accept nodes. */
	std::vector<Node_ptr> getAcceptNodes() const {
		return myacceptNodes;
	}

	/** Return a vector of initial nodes. */
	std::vector<Node_ptr> getInitialNodes() const {
		return myinitialNodes;
	}

	/** Return the vector of Nodes in an automaton. */
	std::vector<Node_ptr> getNodeList() const {
		return myNodeList;
	}

	/** Return node id's. */
	std::vector<new_node_id::node_id_type> getNodeIds() const {
		return myNodeIds;
	}

	/** Return label id's from the alphabet. */
	std::set<label_cache::label_id_type> getAlphabetIds() const {
		return myalphabetIds;
	}

	/** Set labels id's. */
	void setAlphabetIds(const std::set<label_cache::label_id_type> labelIds) {
		myalphabetIds.clear();
		myalphabetIds = labelIds;
	}

};
}
#include "Automaton.hpp"
#endif /* AUTOMATON_H_ */
