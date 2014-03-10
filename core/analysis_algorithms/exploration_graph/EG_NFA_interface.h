/*
 * EG_NFA_interface.h
 *
 *  Created on: Jan 13, 2012
 *      Author: goyal
 */

#ifndef EG_NFA_INTERFACE_H_
#define EG_NFA_INTERFACE_H_

#include "ExGraph.h"
#include "ExGraphTransition.h"
#include "extern/finite_automaton/DFAutomaton.h"
using namespace finite_automaton;

namespace exploration_graph {

/** forward declarations */
typedef ExGraph::shortest_path_info_type shortest_path_info_type;

/** Container type for ExGraphTransition to NFA Labels map. */
typedef std::map<ExGraphTransition_Info_ptr, label_id_type>
		ExTransition_to_label_container_type;

/** Type of abstraction function map local to each exploration graph. */
typedef std::map<ExNode_ptr, ExTransition_to_label_container_type>
		local_abs_function_map_type;

/** Global map type that maps an ExGraph to its local_abs_function. */
typedef std::map<ExGraph*, local_abs_function_map_type> global_map_type;

/** Interface class for NFA and ExGraph. */
class EG_NFA_interface {

public:

	/** Convert exGraph to an NFA given EG_path_info in an exploration graph.
	 *
	 * EG_path_info consists of an initial node pointer and a list of transitions
	 * on that path.
	 *
	 * @note The path is assumed to be without loops.
	 *
	 * All DISCRETE transitions in an EG_path are labeled as "1"; whereas,
	 * rest of the DISCRETE transitions' labeling starts from 2 in a non-decreasing
	 * order. TIME_ELAPSE transition in an EG_path are labeled as "time_bad", while
	 * the rest are labeled as "time_not_bad", and CONTAINMENT as "EPS".
	 *
	 * O((n+m)log n + m log m).
	 */
	static NFAutomaton to_language(const ExGraph& exGraph,
			shortest_path_info_type EG_path_info);

	/** Create DFA from a list of transitions. */
	static DFAutomaton create_DFA_from_path(
			const std::vector<ExGraphTransition> myTransitions);

	/** Compose exGraph with DFA and return the resulting EG. */
	static ExGraph
			from_language(const ExGraph& exGraph, const DFAutomaton& DFA);

	/** Compose exGraph with NFA and return the resulting EG. */
	static ExGraph
			from_language(const ExGraph& exGraph, const NFAutomaton& NFA);

	static global_map_type global_map; //!< global map

	/** Initialize global_map. */
	static global_map_type init_global_map();
};
}
#include "EG_NFA_interface.hpp"
#endif /* EG_NFA_INTERFACE_H_ */
