/*
 * composition_operators.cpp
 *
 *  Created on: Sep 3, 2009
 *      Author: frehse
 */

#include "core/hybrid_automata/composition_operators.h"

#include "core/continuous/continuous_dynamics/continuous_dynamics_composition.h"
//#include "../abstract_framework/continuous/continuous_set_transforms/continuous_set_transform.h""
#include "core/continuous/continuous_set_transforms/continuous_set_transform_composition.h"
#include "core/symbolic_states/symbolic_state.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/hybrid_automata/hybrid_automaton.h"

/** Forward declarations */
namespace continuous {
class continuous_set;
typedef boost::shared_ptr<const continuous_set> continuous_set_const_ptr;
continuous_set_ptr compute_intersection(const continuous_set_const_ptr& p1,
		const continuous_set_const_ptr& p2);
}

namespace hybrid_automata {

variable_id_set compute_input_variables(const hybrid_automaton_const_ptr& h1,
		const hybrid_automaton_const_ptr& h2) {
	variable_id_set var1, var2, inp1, inp2, contr1, contr2, new_inp;

	var1 = h1->get_variable_ids();
	var2 = h2->get_variable_ids();
	inp1 = h1->get_input_variables();
	inp2 = h2->get_input_variables();
	contr1 = var1;
	set_difference_assign(contr1, inp1);
	contr2 = var2;
	set_difference_assign(contr2, inp2);

	// We suppose var_Oi=var_Ci=contri
	// var_O=var_O1 \cup var_O2
	// var_I=(var_I1 \cup var_I2) \ var_O
	new_inp = inp1;
	new_inp.insert(inp2.begin(), inp2.end());
	set_difference_assign(new_inp, contr1);
	set_difference_assign(new_inp, contr2);

	return new_inp;
}

variable_id_set compute_const_variables(const hybrid_automaton_const_ptr& h1,
		const hybrid_automaton_const_ptr& h2) {
	variable_id_set var1, var2, const1, const2, new_const;

	var1 = h1->get_variable_ids();
	var2 = h2->get_variable_ids();
	const1 = h1->get_const_variables();
	const2 = h2->get_const_variables();

	new_const = const1;
	new_const.insert(const2.begin(), const2.end());

	return new_const;
}

location::ptr compose_location(const location::const_ptr& p1, const location::const_ptr& p2) {
	std::string new_name = p1->get_name() + "~" + p2->get_name();

	const time_constraints& tc1 = p1->get_time_constraints();
	const time_constraints& tc2 = p2->get_time_constraints();

	assert(tc1.get_invariant());
	assert(tc2.get_invariant());
	continuous::continuous_set::const_ptr new_inv = continuous::compute_intersection(
			tc1.get_invariant(), tc2.get_invariant());
	if (!new_inv) {
		std::cerr << "inv1:" << tc1.get_invariant() << std::endl;
		std::cerr << "inv2:" << tc2.get_invariant() << std::endl;
		throw std::runtime_error("could not parallel compose invariants!");
	}
	continuous::continuous_dynamics::const_ptr dyn = continuous::parallel_compose(
			tc1.get_dynamics(), tc2.get_dynamics());
	if (!dyn) {
		throw std::runtime_error("could not parallel compose dynamics!");
	}

	time_constraints new_tcons(new_inv, dyn);

	location::ptr new_loc = location::ptr(new location(new_name, new_tcons));
	return new_loc;
}

jump_constraints compose_jump_constraints(const jump_constraints& c1, const jump_constraints& c2) {
	// intersect guard
	continuous::continuous_set::const_ptr g;
	const continuous::continuous_set::const_ptr& g1 = c1.get_guard();
	const continuous::continuous_set::const_ptr& g2 = c2.get_guard();

	if (g1) {
		if (g2) {
			g = continuous::compute_intersection(g1, g2);
		} else {
			g = g1;
		}
	} else {
		g = g2;
	}

	// compose transform
	continuous::continuous_set_transform_ptr t = continuous::parallel_compose(c1.get_transform(),
			c2.get_transform());

	if (!t) {
		throw std::runtime_error("could not parallel compose transitions");
	}

	return jump_constraints(g, t);
}

jump_constraints compose_autonomous_jump_constraints(const jump_constraints& c1,
		const hybrid_automaton_const_ptr& h1, const hybrid_automaton_const_ptr& h2) {
	variable_id_set contr_vars1 = h1->get_controlled_variables();
	variable_id_set contr_vars2 = h2->get_controlled_variables();


	// For non-I/O semantics, don't force constness on variables that are controlled by h1 as well.
	// In an I/O framework, variables are controlled by exactly one automaton, so
	// this rule doesn't interfere with I/O.
	variable_id_set const_vars = contr_vars2;
	set_difference_assign(const_vars, contr_vars1);

	jump_constraints jcons = c1;
	if (!const_vars.empty()) {
		continuous::continuous_set_transform_ptr t = extend_transform_with_const_variables(
				jcons.get_transform(), const_vars);
		//assert(t);
		jcons.set_transform(t);
	}
	return jcons;
}

symbolic_state_collection_ptr compose_symbolic_states(
		const symbolic_state_collection_const_ptr& s1,
		const symbolic_state_collection_const_ptr& s2) {
	// pairwise enumerate the states and intersect the symbolic state
	// @todo this should just be the intersection of the symbolic_state_collections

	symbolic_state_collection_ptr res;
	if (!s1) {
		if (!s2) {
			// case !s1 !s2
			res = symbolic_state_collection_ptr();
		} else {
			// case !s1 s2
			res = symbolic_state_collection_ptr(s2->clone());
		}
	} else {
		if (!s2) {
			// case s1 !s2
			res = symbolic_state_collection_ptr(s1->clone());
		} else {
			// case s1 s2
			// Create new set
			res=symbolic_state_collection_ptr(s1->create());
			for (symbolic_state_collection::const_iterator i1 = s1->begin(); i1
					!= s1->end(); ++i1) {
				for (symbolic_state_collection::const_iterator i2 = s2->begin(); i2
						!= s2->end(); ++i2) {
					symbolic_state::ptr s = *i1;
					s->intersection_assign(*i2);
					if (!math::definitely(s->is_empty()))
						res->add(s);
				}
			}
		}
	}
	return res;
}

}


