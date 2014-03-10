#ifndef NFAUTOMATON_H_
#define NFAUTOMATON_H_

/*********************************************************************
 * NFAutomaton.h													 *
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

#include "Automaton.h"

namespace exploration_graph {
class EG_NFA_interface;
}

namespace finite_automaton {

/** forward declarations */
class DFAutomaton;

class NFAutomaton: public Automaton, public boost::enable_shared_from_this<
		NFAutomaton> {

public:

	typedef boost::shared_ptr<NFAutomaton> NFA_ptr;

	/** Return a shared_ptr to *this. */
	NFA_ptr get_ptr() {
		return boost::enable_shared_from_this<NFAutomaton>::shared_from_this();
	}
	;

	/** NFAutomaton Constructor */
	NFAutomaton();

	/** Construct a deep copy of automaton. */
	NFAutomaton(const NFAutomaton& automaton);

	/** NFAutomaton destructor */
	virtual ~NFAutomaton();

	/** Convert an NFA into a DFA.
	 * O(2^n(n log n)).*/
	DFAutomaton convertToDFA() const;

	/** Compute Kleene's star of a NFA.
	 * O(n+m)log n.*/
	NFAutomaton computeKleeneStar() const;

	/** Concatenate an NFA with secondNFA.
	 * O((n1 + m1)log n1 + (n2 + m2)log n2).*/
	NFAutomaton concatenate(const NFAutomaton& secondNFA) const;

	/** Bellman ford algorithm is used.
	 * Multiple sources, if any, are merged first and then
	 * separated. */
	virtual std::vector<std::string> getShortestPath() const;

private:

	/** flag = true meant for operations such as concatenate,
	 * KleeneStar and convertToDFA.
	 * flag = false is meant for getShortestPath(). */
	void mergeInitialNodes(bool flag) {

		std::vector<Node_ptr> initialNodes = this->getInitialNodes();
		std::set<label_cache::label_id_type> labelIds = this->getAlphabetIds();

		if (initialNodes.size() > 1) {

			label_cache::label_id_type lId;
			if (flag == false)
				lId = *(labelIds.begin());
			else if (flag == true) {
				lId = label_cache::label_id_type(0);
				this->addToAlphabet("EPS");
			}
			const AutomatonNode& newNode = this->createNode();
			Node_ptr newNode_shared_pointer = access_method::access(newNode);
			for (std::vector<int>::size_type i = 0; i < initialNodes.size();
					++i) {
				newNode_shared_pointer->addSuccNode(*(initialNodes[i]), 1, lId);
			}
			this->resetInitialNode(newNode_shared_pointer);
		} //!< Now, we have just one NFA initial state
		return;
	}

	/** Separate initial nodes which were merged before. */
	void separateInitialNodes() {

		std::vector<Node_ptr> initialNodes = this->getInitialNodes();
		Node_ptr newInitialNode = initialNodes.at(0);

		initialNodes.clear();
		std::vector<AutomatonEdge> myEdges = newInitialNode->getEdges();

		for (std::vector<AutomatonEdge>::iterator eIt = myEdges.begin();
				eIt != myEdges.end(); ++eIt) {

			const AutomatonNode& succNodeObj_ref = (*eIt).getSuccessorNode();
			Node_ptr succNode = access_method::access(succNodeObj_ref);
			initialNodes.push_back(succNode);
		}

		std::vector<Node_ptr>::iterator initNodes_it;
		initNodes_it = initialNodes.begin();
		this->resetInitialNode(*initNodes_it);
		initNodes_it++;
		while (initNodes_it != initialNodes.end()) {
			this->setInitialNode(*(*initNodes_it));
			initNodes_it++;
		}

		std::vector<Node_ptr> NodeList = this->getNodeList();
		std::vector<new_node_id::node_id_type> NodeIds = this->getNodeIds();

		NodeList.erase(NodeList.begin() + (NodeList.size() - 1));
		NodeIds.erase(NodeIds.begin() + (NodeIds.size() - 1));
		this->resetNodeList(NodeList);
		this->resetNodeIds(NodeIds);

		return;
	}

	friend class DFAutomaton;
	friend class exploration_graph::EG_NFA_interface;
};
}

#include "DFAutomaton.h"
#include "NFAutomaton.hpp"
#endif /* NFAUTOMATON_H_ */
