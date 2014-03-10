#ifndef REACHABILITY_ALGORITHM_H_
#define REACHABILITY_ALGORITHM_H_

#include "boost/shared_ptr.hpp"

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
typedef boost::shared_ptr<const hybrid_automaton> hybrid_automaton_const_ptr;
class symbolic_state;
typedef boost::shared_ptr<symbolic_state> symbolic_state_ptr;
class symbolic_state_collection;
typedef boost::shared_ptr<symbolic_state_collection> symbolic_state_collection_ptr;
typedef boost::shared_ptr<const symbolic_state_collection> symbolic_state_collection_const_ptr;
class post_operator;
typedef boost::shared_ptr<post_operator> post_operator_ptr;
class discrete_post;
typedef boost::shared_ptr<discrete_post> discrete_post_ptr;
typedef boost::shared_ptr<const discrete_post> discrete_post_const_ptr;
class continuous_post;
typedef boost::shared_ptr<continuous_post> continuous_post_ptr;
typedef boost::shared_ptr<const continuous_post>
		continuous_post_const_ptr;
class passed_and_waiting_list;
typedef boost::shared_ptr<passed_and_waiting_list> passed_and_waiting_list_ptr;
class reachability_scenario;
}

namespace hybrid_automata {

/**
 * The base class providing the reachability algorithm for hybrid automaton. Implementation of this class
 * are supposed to override the reachability algorithm accordingly.
 */

class reachability_algorithm {
public:
	/*! \brief
	 * C++ class Constructor.
	 */
	reachability_algorithm(const reachability_scenario& scen);
	/*! \brief
	 * C++ class destructor.
	 */
	virtual ~reachability_algorithm();

	/** Returns the collection of symbolic states which are reachable in the hybrid automaton h starting
	 * from start_states. Different variations of the reachability can be implemented by overriding
	 * the base implementation in the derived classes. The variations can be in the order in which the post
	 * operators are applied, depth first or breadth first traversal, exact computation or overapproximation
	 * etc.
	 */
	virtual symbolic_state_collection_ptr reach(hybrid_automaton_ptr H,
			const symbolic_state_collection_const_ptr& start_states);

	/** Returns the collection of symbolic states which are reachable in the hybrid automaton h starting
	 * from the initial states.
	 */
	virtual symbolic_state_collection_ptr reach(hybrid_automaton_ptr H);

protected:
    void apply_post_ops(const symbolic_state_ptr & psstate);

    /* Returns true if next iteration should be carried out and false
     * otherwise. */
    bool iter_check(int iter_count);

	const reachability_scenario& my_scen;
	passed_and_waiting_list_ptr my_plwl;
	hybrid_automaton_ptr my_aut;
	continuous_post_ptr my_post_c;
	discrete_post_ptr my_post_d;
};

}

#endif /*REACHABILITY_ALGORITHM_H_*/
