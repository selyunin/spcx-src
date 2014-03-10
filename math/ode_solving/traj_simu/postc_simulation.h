/*
 * continuous_post.h
 *
 *  Created on: Jun 23, 2009
 *      Author: frehse
 */

#ifndef CONTINUOUS_POST_SIMULATION_H_
#define CONTINUOUS_POST_SIMULATION_H_

#include "core/post_operators/post_operator.h"
#include "core/post_operators/continuous_post.h"
#include "core/continuous/polyhedra/hyperbox/finite_hyperbox.h"
#include "core/symbolic_states/symbolic_state_collection.h"
#include "core/symbolic_states/symbolic_state.h"
#include "core/symbolic_states/symbolic_state_collection_stl_list.h"

#include "math/ode_solving/ode_solver.h"

#include "continuous_set_simulation.h"

#include <boost/random.hpp>

/** Forward declaration of classes used in header file. */

namespace hybrid_automata {

/** Continuous post operator for the simulation scenario
 *
 * When applied to the initial set of states,
 * the operator samples the set to produce a set of
 * states for processing by the discrete post operator.
 *
 * Afterwards, in the fixpoint loop, the operator has no function:
 * it returns the set (a set of points) unchanged.
 */

class postc_simulation: public continuous_post {
public:
	typedef boost::shared_ptr<postc_simulation> ptr;
	typedef boost::shared_ptr<const postc_simulation> const_ptr;

	//TODO
	typedef double scalar_type;
	typedef math::vector<scalar_type> vector;
	typedef math::vdom_vector<scalar_type> state;
	typedef continuous::continuous_set_simulation<scalar_type>::trajectory trajectory;

	typedef math::ode::ode_solver<scalar_type>::ode_parameters ode_parameters;

	postc_simulation(){
		my_time_horizon = -1;
		my_params = ode_parameters();
		my_local_time_horizon = -1;
		my_bbox= continuous::finite_hyperbox<scalar_type>();
		my_uniform_sampling = 0;
		my_rng = boost::mt19937();
	}


	postc_simulation(const postc_simulation & p) {
		my_time_horizon = p.get_global_time_horizon();
		my_params = p.get_ode_parameters() ;
		my_local_time_horizon = p.get_time_horizon();
		my_bbox=p.get_bounding_box();
		my_uniform_sampling = p.get_init_uniform_sampling();
		my_rng = boost::mt19937();
	};

	virtual ~postc_simulation() {
	}
	;
	virtual continuous_post* clone(){
		return new postc_simulation(*this);

	}

	/** Return the continuous_set that results from applying time elapse
	 * with the time constraints tcons to the continuous set *cset.
	 *
	 * If the returned set is empty, the implementation must return
	 * a null pointer (continuous_set_ptr()). This is to avoid
	 * emptiness checks further down the line.
	 *
	 * This is the only function (apart from clone()) that needs to be provided by derived
	 * implementations.*/
	virtual continuous::continuous_set_ptr post(const time_constraints& tcons,
			const continuous::continuous_set_const_ptr& cset) const{

		throw std::runtime_error("feature not implemented yet");

	}

	/** Return the continuous_set that results from applying time elapse
	 * in the discrete set *dset (all locations of which are considered
	 * equivalent w.r.t. time elapse) to the continuous set *cset.
	 *
	 * The default implementation simply chooses the time constraints of
	 * the first location in the set.*/
	virtual continuous::continuous_set_ptr post(
			const hybrid_automaton_const_ptr& aut,
			const discrete::discrete_set_const_ptr& dset,
			const continuous::continuous_set_const_ptr& cset) const{

		throw std::runtime_error("feature not implemented yet");


	}

	using post_operator::add_post_states;

	/** This continuous post operator just passes the state to the waiting result set. */
	virtual void add_post_states(const hybrid_automaton_const_ptr& aut,
			symbolic_state_collection_ptr& passed_result_set,
			symbolic_state_collection_ptr& waiting_result_set,
			const symbolic_state_ptr& sstate) const;

	/** This continuous post operator just passes the states to the waiting result set. */
	virtual void add_post_states(const hybrid_automaton_const_ptr& aut,
			symbolic_state_collection_ptr& passed_result_set,
			symbolic_state_collection_ptr& waiting_result_set,
			const symbolic_state_collection_const_ptr& sstate_set) const;

	virtual symbolic_state_collection_ptr select_states_from_polyhedron(const symbolic_state_ptr& poly) const;

	virtual state pick_uniform_state_from_hyperbox(const continuous::hyperbox<scalar_type>& hyp) const;

	virtual scalar_type gen_random_scalar(scalar_type min, scalar_type max) const;

	/** Set time horizon to th.
	 *
	 * A negative value denotes unlimited time. */
	virtual void set_global_time_horizon(double th){ my_time_horizon = th;}
	virtual double get_global_time_horizon() const { return my_time_horizon; }

	 /* A negative value denotes unlimited time. */
	virtual void set_time_horizon(double lth){my_local_time_horizon = lth;}
	virtual double get_time_horizon() const  { return my_local_time_horizon; }


	/** Set sampling time to ts.
	 *
	 * A negative value denotes no sampling. */
	virtual void set_ode_parameters(ode_parameters ts) {my_params = ts;}
	virtual ode_parameters get_ode_parameters() const  { return my_params; }


	virtual void set_bounding_box(continuous::finite_hyperbox<scalar_type> hbox){
		my_bbox = hbox;
	}
	virtual const continuous::finite_hyperbox<scalar_type> & get_bounding_box() const {
			return my_bbox;
	}

	virtual void set_init_uniform_sampling(unsigned int s){
		my_uniform_sampling = s;
	}
	virtual unsigned int get_init_uniform_sampling() const {
			return my_uniform_sampling;
		}


private:
	double my_time_horizon;
	double my_local_time_horizon;
	unsigned int my_uniform_sampling;
	ode_parameters my_params;
	continuous::finite_hyperbox<scalar_type> my_bbox;
	mutable boost::mt19937 my_rng;
};

}

#endif /* CONTINUOUS_POST_H_ */
