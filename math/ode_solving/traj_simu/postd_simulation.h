

#ifndef DISCRETE_POST_SIMULATION_H_
#define DISCRETE_POST_SIMULATION_H_

#include "core/post_operators/discrete_post.h"
#include "math/ode_solving/traj_simu/postc_simulation.h"

#include "core/continuous/continuous_set_collection.h"
#include "core/post_operators/post_operator.h"

#include "core/hybrid_automata/named_label.h"
#include "core/hybrid_automata/transition_id.h"

#include "core/continuous/continuous_dynamics/continuous_dynamics.h"
#include "core/pwl/passed_and_waiting_list.h"
#include "math/ode_solving/traj_simu/traj_simu.h"

namespace hybrid_automata {

/** A discrete post operator is a post operator that iterates
 * over a set of transitions provided by the automaton.
 *
 * An implementation must provide two methods:
 * - clone() creates a copy of the operator
 * - post(jmp,cset) computes the continuous set resulting from applying
 *   the jump constraints jmp to the continuous set cset.
 */

class postd_simulation: public discrete_post {
public:
	typedef boost::shared_ptr<postd_simulation> ptr;
	typedef boost::shared_ptr<const postd_simulation> const_ptr;

	typedef double scalar_type;

	typedef traj_simu<scalar_type> simu;

	postd_simulation(const postc_simulation::const_ptr & postc) :
		mypostc(*postc) {

		mypostcptr = postc;
		simulation = traj_simu<scalar_type>(mypostc.get_ode_parameters(),
							traj_simu<scalar_type>::traj_simu_parameters(
									mypostc.get_global_time_horizon(),mypostc.get_time_horizon()) );
	}

	postd_simulation(const postd_simulation & postd) :
		mypostc(*postd.get_postc()) {
		mypostcptr = postd.get_postc();
		simulation = traj_simu<scalar_type>(mypostc.get_ode_parameters(),
									traj_simu<scalar_type>::traj_simu_parameters(
											mypostc.get_global_time_horizon(),mypostc.get_time_horizon()) );

		my_simu_algo = postd.my_simu_algo;
	}

	virtual ~postd_simulation() {
	}
	;

	void reset_options(){

		simulation = traj_simu<scalar_type>(mypostc.get_ode_parameters(),
									traj_simu<scalar_type>::traj_simu_parameters(
											mypostc.get_global_time_horizon(),mypostc.get_time_horizon()) );

	}

	void set_simu_internal_algo(unsigned int algo){
		this->my_simu_algo = algo;
	}

	unsigned int get_simu_internal_algo() const {
		return my_simu_algo;
	}

	postc_simulation::const_ptr get_postc() const { return mypostcptr; }

	/** Shallow copy. */
	virtual discrete_post* clone() ;

	/** Return the continuous_set that results from applying the transition
	 * *trans to the continuous set *cset.
	 *
	 * This is the only function that needs to be provided by derived
	 * implementations.*/
	virtual continuous::continuous_set_collection post(
			const jump_constraints& trans,
			continuous::continuous_set_const_ptr source_inv,
			continuous::continuous_set_const_ptr target_inv,
			continuous::continuous_set_const_ptr cset) const;//TODO pure virt

	/** Apply the post operator for the transitions with trans_id to the symbolic
	 * state sstate and add the result to result_set.
	 *
	 * The target discrete set is the target location of trans. The
	 * resulting continuous set or sets are obtained by calling post(trans,target_inv,cset).
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

	virtual void add_post_states(const hybrid_automaton_const_ptr& aut,
			symbolic_state_collection_ptr& passed_result_set,
			symbolic_state_collection_ptr& waiting_result_set,
			const std::list<transition_id>& trans_set, const symbolic_state_ptr& sstate) const;

	virtual std::vector<continuous::continuous_set::ptr>  post(
			const std::vector<jump_constraints> & trans,
			const continuous::continuous_set_const_ptr & source_inv,
			const continuous::continuous_dynamics::const_ptr & source_dyn,
			const std::vector<continuous::continuous_set_const_ptr> & targets_inv,
			const continuous::continuous_set_const_ptr & cset,
			continuous::continuous_set_simulation<scalar_type>::ptr& cset_trajs) const;


private:
	postc_simulation::const_ptr mypostcptr;
	const postc_simulation & mypostc;
	traj_simu<scalar_type> simulation;
	unsigned int my_simu_algo;
};

}

#endif /* DISCRETE_POST_H_ */
