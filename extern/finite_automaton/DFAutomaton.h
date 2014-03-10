#ifndef DFAUTOMATON_H_
#define DFAUTOMATON_H_

/*********************************************************************
 * DFAutomaton.h													 *
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
class NFAutomaton;
typedef AutomatonEdge::cost_type distance_type;

/**
 * DFA class represents a Deterministic Finite Automaton.
 *
 * For consistency, we deal with complete DFAs. A complete
 * DFA has exactly one outgoing transition from each node
 * for every label in its alphabet. Thus, a DFA is completed
 * first if not already.
 *
 * Furthermore, before carrying out a binary operation on two
 * DFAs, alphabet of each DFA is extended with respect to the
 * alphabet of other DFA and then, each DFA is completed with
 * its new alphabet. The complete process is achieved with
 * extend() method.
 */

class DFAutomaton: public Automaton, public boost::enable_shared_from_this<
		DFAutomaton> {

public:

	typedef boost::shared_ptr<DFAutomaton> DFA_ptr;

	/** Return a shared_ptr to *this. */
	DFA_ptr get_ptr() {
		return boost::enable_shared_from_this<DFAutomaton>::shared_from_this();
	}
	;

	/** DFAutomaton Constructors */
	DFAutomaton();

	/** Create a deep copy of automaton. */
	DFAutomaton(const DFAutomaton& automaton);

	/** DFAutomaton destructor */
	~DFAutomaton();

	/** Return the shortest word in a DFA.
	 *
	 * @note Its implemented using BFS, which is an optimal solution when
	 * all the edge costs are same. O(n+m)log n.
	 *
	 * @cite http://en.wikipedia.org/wiki/Breadth-first_search
	 */
	std::vector<std::string> getShortestWord() const;

	/**
	 * Minimize a DFA
	 *
	 * Hopcroft's partition refinement approach is used.
	 * @cite http://en.wikipedia.org/wiki/DFA_minimization
	 * O(l(n*n)log n).
	 */
	DFAutomaton minimize() const;

	/**
	 * The method converts an incomplete DFA into a complete DFA.
	 *
	 * A new node is added which has transitions to itself on each
	 * label in its alphabet. Also, all the undefined transitions
	 * from other nodes are mapped to newly added node.
	 * O(l + m log l).*/
	void makeComplete();

	/**
	 * Compute and return the complement of a DFA.
	 * @note Accepting and non-accepting nodes are interchanged.
	 * O(n+m)log n.
	 */
	DFAutomaton computeComplement() const;

	/** Compute and return the Union of a DFA with secondDFA.
	 * O((n1 * n2)(l1 + l2) log (n1 * n2)).*/
	DFAutomaton computeUnion(const DFAutomaton& secondDFA) const;

	/** Compute and return the Intersection of a DFA with secondDFA.
	 * O((n1 * n2)(l1 + l2) log (n1 * n2)).*/
	DFAutomaton computeIntersection(const DFAutomaton& secondDFA) const;

	/** Compute and return the Union of a DFA with secondDFA.
	 * O((n1 * n2)(l1 + l2) log (n1 * n2)).*/
	DFAutomaton computeDifference(const DFAutomaton& secondDFA) const;

	/** Concatenate a DFA with secondDFA.
	 * O(2^(n1+n2)((n1+n2) log (n1+n2))).*/
	DFAutomaton concatenate(const DFAutomaton& secondDFA) const;

	/** Convert a DFA into an NFA.
	 *
	 * The method copies the whole DFA structure into its NFA counterpart.
	 * O(n+m)log n. */
	NFAutomaton convertToNFA() const;

	/**	Compute Kleene's star of a DFA.
	 * 	@note Kleene's star is an NFA.	O(n+m)log n.*/
	NFAutomaton computeKleeneStar() const;

	/** Return true if the word is recognized by this DFA language
	 * else, return false. O(n+m).*/
	bool Accepts(const std::string word) const;

protected:

	/** First, extend the alphabet of a DFA with respect to
	 * the given alphabet and then, completes the DFA. */
	void extend(const std::set<label_cache::label_id_type> LabelIds);

	/** Convert a DFA to NFAStar.
	 *
	 * Self EPS transitions are added to each node.
	 * Used during composition of an exploration graph
	 * with an NFA. */
	NFAutomaton convertToNFAStar() const;

	friend class exploration_graph::EG_NFA_interface;
	friend class NFAutomaton;

};
}
#include "NFAutomaton.h"
#include "DFAutomaton.hpp"

#endif /* DFAUTOMATON_H_ */
