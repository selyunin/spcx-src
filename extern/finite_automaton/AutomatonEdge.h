#ifndef AUTOMATONEDGE_H_
#define AUTOMATONEDGE_H_

/*********************************************************************
 * AutomatonEdge.h													 *
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

#include "label_cache.h"

namespace exploration_graph {
class EG_NFA_interface;
}

namespace finite_automaton {

class AutomatonNode;
typedef boost::shared_ptr<AutomatonNode> Node_ptr;

/**
 * AutomatonEdge class represents the outgoing transition
 * from an automaton node. A transition is defined by a
 * pointer to the successor node, the cost and label id.
 */
class AutomatonEdge {
public:

	typedef label_cache::label_id_type label_id_type;

	typedef unsigned int cost_type;

	/** AutomatonEdge destructor */
	~AutomatonEdge();

	/** AutomatonEdge constructors */
	AutomatonEdge();

	/** Return reference to the successor node */
	const AutomatonNode& getSuccessorNode() const;

	/** Return edge cost */
	const cost_type& getEdgeCost() const;

	/** Return edge label */
	const std::string getEdgeLabel() const;

protected:

	/** Add successor node */
	AutomatonEdge(const AutomatonNode& eSuccNode, const cost_type& eCost,
			const label_id_type& eLabelId);

private:

	label_id_type edgeLabel_id; /*!< edge label id */

	Node_ptr successorNode; /*!< pointer to the successor node */

	cost_type edgeCost; /*!< cost of the edge */

	/** Return edge label id */
	const label_id_type& getEdgeLabelId() const {
		return edgeLabel_id;
	}

	/** These classes make use of getEdgeLabelId() method. */
	friend class Automaton;
	friend class DFAutomaton;
	friend class NFAutomaton;
	friend class AutomatonNode;

	friend class exploration_graph::EG_NFA_interface;

};
}
#include "AutomatonNode.h"
#include "AutomatonEdge.hpp"
#endif /* AUTOMATONEDGE_H_ */
