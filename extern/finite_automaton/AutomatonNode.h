#ifndef AUTOMATONNODE_H_
#define AUTOMATONNODE_H_

/*********************************************************************
 * AutomatonNode.h													 *
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

#include "boost/shared_ptr.hpp"

namespace finite_automaton {

/** forward declarations */
class AutomatonEdge;
typedef AutomatonEdge::cost_type cost_type;

/**
 * AutomatonNode class represents the node of an automaton.
 * Each node is defined by a vector of edges, which denotes
 * the outgoing transitions from this node.
 */
class AutomatonNode: public boost::enable_shared_from_this<AutomatonNode> {
public:

	typedef boost::shared_ptr<AutomatonNode> Node_ptr;

	/** Return a shared_ptr to *this. */
	Node_ptr get_ptr() {
		return boost::enable_shared_from_this<AutomatonNode>::shared_from_this();
	}
	;

	/** Node status type*/
	typedef enum {
		IS_VISITED, IS_NOT_VISITED
	} status;

	/** AutomatonNode destructor */
	~AutomatonNode();

	/** AutomatonNode constructor */
	AutomatonNode();

	/** Return the vector of outgoing transitions for this node. */
	std::vector<AutomatonEdge> getEdges() const;

	/** Get node's status. */
	status getStatus() const;

	/** Set node's status. */
	void setStatus(const status st);

private:

	/** Add a new edge.
	 *
	 * @note label is stored as an id while constructing the edge */
	void addSuccNode(const AutomatonNode& succNode, const cost_type& cost,
			const label_cache::label_id_type& label_id) {
		AutomatonEdge Edge(succNode, cost, label_id);
		myEdges.push_back(Edge);
	}

	/** Set the successor node at index Idx in the edge vector.
	 *
	 * @note Used during DFA minimization. */
	void setSuccNodeAtIdx(const AutomatonNode& succNode, int Idx) {
		AutomatonEdge Edge(succNode, myEdges[Idx].getEdgeCost(),
				myEdges[Idx].getEdgeLabelId());
		myEdges[Idx] = Edge;
	}

	std::vector<AutomatonEdge> myEdges; /*!< list of outgoing transitions of the node */

	status myStatus; /*!< Node status. */

	/** It uses addSuccNode() method. */
	friend class Automaton;

	/** It uses addSuccNode() and setSuccNodeAtIdx() methods. */
	friend class DFAutomaton;

	/** It uses addSuccNode() method. */
	friend class NFAutomaton;
};

}
#include "AutomatonEdge.h"
#include "AutomatonNode.hpp"

#endif /* AUTOMATONNODE_H_ */
