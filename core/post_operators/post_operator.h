#ifndef POST_OPERATOR_H_
#define POST_OPERATOR_H_

#include "boost/shared_ptr.hpp"
//#include "../hybrid_automata/hybrid_automaton.h"

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
typedef boost::shared_ptr<const hybrid_automaton> hybrid_automaton_const_ptr;
class symbolic_state;
typedef boost::shared_ptr<symbolic_state> symbolic_state_ptr;
typedef boost::shared_ptr<const symbolic_state> symbolic_state_const_ptr;
class symbolic_state_collection;
typedef boost::shared_ptr<symbolic_state_collection> symbolic_state_collection_ptr;
typedef boost::shared_ptr<const symbolic_state_collection> symbolic_state_collection_const_ptr;
}

namespace hybrid_automata {

/** A post operator is an operator that adds the result of an operation on
 * a symbolic states sstate to a symbolic state collection result_set.
 */

class post_operator {
public:
	typedef boost::shared_ptr<post_operator> ptr;
	typedef boost::shared_ptr<const post_operator> const_ptr;

	virtual ~post_operator() {
	}
	;

	/** Virtual copy constructor. */
	virtual post_operator* clone() = 0;

	/** Apply the post operator to the symbolic state sstate and
	 * add the result to result_set. */
	virtual void add_post_states(const hybrid_automaton_const_ptr& aut,
			symbolic_state_collection_ptr& passed_result_set,
			symbolic_state_collection_ptr& waiting_result_set,
			const symbolic_state_ptr& sstate) const = 0;

	/** Apply the post operator to all symbolic states in sstate_set and
	 * add the result to result_set.
	 *
	 * The default implementation enumerates the symbolic states in sstate_set
	 * and calls add_post_states for each one.
	 */
	virtual void add_post_states(const hybrid_automaton_const_ptr& aut,
			symbolic_state_collection_ptr& passed_result_set,
			symbolic_state_collection_ptr& waiting_result_set,
			const symbolic_state_collection_const_ptr& sstate_set) const;
};

}

#endif /*POST_OPERATOR_H_*/
