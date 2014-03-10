#include "core/post_operators/post_operator.h"
#include "core/symbolic_states/symbolic_state_collection.h"

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
typedef boost::shared_ptr<const hybrid_automaton> hybrid_automaton_const_ptr;
}

namespace hybrid_automata {

void post_operator::add_post_states(const hybrid_automaton_const_ptr& aut,
		symbolic_state_collection_ptr& passed_result_set,
		symbolic_state_collection_ptr& waiting_result_set,
		const symbolic_state_collection::const_ptr& sstate_set) const {
	for (symbolic_state_collection::const_iterator it=
			sstate_set->begin(); it !=sstate_set->end(); ++it) {
		add_post_states(aut, passed_result_set, waiting_result_set, *it);
	}
}

}
