/*
 * convert_predicate.cpp
 *
 *  Created on: Sep 7, 2009
 *      Author: frehse
 */

#include "core/continuous/convert_predicate.h"

//#include "../abstract_framework/hybrid_automata/location_constraint_set.h"
#include "core/discrete/singleton_set.h"
#include "core/hybrid_automata/location_eq_node.h"
#include "core/symbolic_states/symbolic_state.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/predicates/share_continuous_and_discrete_tree.h"
#include "core/continuous/predicate_continuous_set_constructors.h"
#include "core/hybrid_automata/location_constraint_set_constructors.h"

//#include "ppl_polyhedron/continuous_set_PPL_NNC.h"
//#include "ppl_polyhedron/convert_to_ppl.h"

namespace hybrid_automata {
void convert_and_add_predicate_to_symbolic_state_collection(const tree::node::ptr& p,
		hybrid_automata::symbolic_state_collection& sstates) {
	// convert the tree to DNF and split into continuous and discrete parts
	typedef std::vector<std::pair<tree::node::ptr, tree::node::ptr> > tree_pair_vector;
	tree_pair_vector res;

	share_continuous_and_discrete_tree(p, res);

	// for each of the pairs, create a symbolic state
	for (tree_pair_vector::iterator it = res.begin(); it != res.end(); ++it) {
		// create the discrete set
		discrete::discrete_set_ptr d = construct_discrete_set_from_conjunction(it->second);

		// create the continuous set
		continuous::continuous_set::ptr c = continuous::construct_predicate_continuous_set(it->first);

		// create the symbolic state
		hybrid_automata::symbolic_state::ptr s = hybrid_automata::symbolic_state::ptr(new hybrid_automata::symbolic_state(d,c));

		sstates.add(s);
	}
}

//continuous_set_ptr convert_affine_predicate_to_continuous_set(
//		const tree::node::ptr& p) {
//	return ppl_polyhedron::convert_to_continuous_set_PPL_NNC(p);
//}

}



