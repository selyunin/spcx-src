/*
 * discrete_post.h
 *
 *  Created on: Jun 23, 2009
 *      Author: frehse
 */

#ifndef DISCRETE_POST_H_
#define DISCRETE_POST_H_

#include "core/continuous/continuous_set_collection.h"
#include "core/post_operators/post_operator.h"

#include "core/hybrid_automata/named_label.h"
#include "core/hybrid_automata/transition_id.h"

/** Forward declaration of classes used in header file. */
namespace continuous {
class continuous_set;
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
typedef boost::shared_ptr<const continuous_set> continuous_set_const_ptr;
}
namespace hybrid_automata {
class jump_constraints;
}

namespace hybrid_automata {

/** A discrete post operator is a post operator that iterates
 * over a set of transitions provided by the automaton.
 *
 * An implementation must provide two methods:
 * - clone() creates a copy of the operator
 * - post(jmp,cset) computes the continuous set resulting from applying
 *   the jump constraints jmp to the continuous set cset.
 */

class discrete_post: public post_operator {
public:
	typedef boost::shared_ptr<discrete_post> ptr;
	typedef boost::shared_ptr<const discrete_post> const_ptr;

	virtual ~discrete_post() {
	}
	;
	virtual discrete_post* clone() = 0;

	/** Return the continuous_set that results from applying the transition
	 * *trans to the continuous set *cset.
	 *
	 * This is the only function that needs to be provided by derived
	 * implementations.*/
	virtual continuous::continuous_set_collection post(
			const jump_constraints& trans,
			continuous::continuous_set_const_ptr source_inv,
			continuous::continuous_set_const_ptr target_inv,
			continuous::continuous_set_const_ptr cset) const = 0;

	/** Pull in declarations from base class */
	using post_operator::add_post_states;

	/** Apply the post operator for the transitions with trans_id to the symbolic
	 * state sstate and add the result to result_set.
	 *
	 * The target discrete set is the target location of trans. The
	 * resulting continuous set or sets are obtained by calling post(trans,target_inv,cset).
	 *
	 * The default discrete post operator sends his states to the waiting_result_set.
	 *
	 * @note post(trans,target_inv,cset) is allowed to return null pointers inside
	 * the collection or a collection with no elements, both signifying an
	 * empty set of states.
	 *
	 * @note It's called add_post_states_trans because the signature of transition_id
	 * is the same as for label_id.
	 * */
	virtual void add_post_states_trans(const hybrid_automaton_const_ptr& aut,
			symbolic_state_collection_ptr& passed_result_set,
			symbolic_state_collection_ptr& waiting_result_set,
			const transition_id& trans, const symbolic_state_ptr& sstate) const;

	/** Apply the post operator over all transitions with label lab to the symbolic
	 * state sstate and add the result to result_set. */
	virtual void add_post_states(const hybrid_automaton_const_ptr& aut,
			symbolic_state_collection_ptr& passed_result_set,
			symbolic_state_collection_ptr& waiting_result_set,
			const label_id& lab,
			const symbolic_state_ptr& sstate) const;

	/** Apply the post operator over the set of labels lab_set to the symbolic
	 * state sstate and add the result to result_set. */
	virtual void add_post_states(const hybrid_automaton_const_ptr& aut,
			symbolic_state_collection_ptr& passed_result_set,
			symbolic_state_collection_ptr& waiting_result_set,
			const label_id_set& lab_set, const symbolic_state_ptr& sstate) const;

	/** Apply the post operator to the symbolic state sstate and
	 * add the result to result_set.
	 *
	 * Iterates over all labels, including the silent label, which is attributed
	 * to transitions without label. */
	virtual void add_post_states(const hybrid_automaton_const_ptr& aut,
			symbolic_state_collection_ptr& passed_result_set,
			symbolic_state_collection_ptr& waiting_result_set,
			const symbolic_state_ptr& sstate) const;
};

}

#endif /* DISCRETE_POST_H_ */
