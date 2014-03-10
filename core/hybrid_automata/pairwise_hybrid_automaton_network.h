/*
 * pairwise_hybrid_automaton_network.h
 *
 *  Created on: Aug 26, 2009
 *      Author: frehse
 */

#ifndef PAIRWISE_HYBRID_AUTOMATON_NETWORK_H_
#define PAIRWISE_HYBRID_AUTOMATON_NETWORK_H_

#include <boost/enable_shared_from_this.hpp>
#include "core/hybrid_automata/hybrid_automaton_network.h"

namespace hybrid_automata {

class pw_network_node;

/** Hybrid automaton implementation to represent a composition of two
 * hybrid automata. It explicitly instantiates the composition as a new
 * hybrid automaton.
 *
 * @note There justification for two classes, pairwise network and hybrid_automaton_pair
 * is that there's not always a pair. Having the network instantiate a pair
 * when necessary (after the first composition) simply avoid awkward tests
 * for null pointers all over. Maybe this is not the optimal solution,
 * since it introduces an additional wrapper, but it seems to look quite clean
 * in almost every other respect.
 */
class pairwise_hybrid_automaton_network: public hybrid_automaton_network {
public:
	typedef boost::shared_ptr<pairwise_hybrid_automaton_network> ptr;
	typedef boost::shared_ptr<const pairwise_hybrid_automaton_network>
			const_ptr;

	// --------------------------------------------
	/** \name Hybrid Automaton Network Interface
	 *  \{ */
	// --------------------------------------------

	/** Create a network consisting of the automaton initial_aut. */
	pairwise_hybrid_automaton_network();

	virtual ~pairwise_hybrid_automaton_network();

	/** Return a shared_ptr to *this. */
	virtual ptr get_ptr();

	/** Return a shared_ptr to const *this. */
	virtual const_ptr get_const_ptr() const;

	/** \brief Return a pointer to a network composed with aut.
	 *
	 * @attention The operation erases the initial states assigned
	 * so far. */
	virtual hybrid_automaton_network::ptr compute_or_assign_composition(
			hybrid_automaton::ptr aut);

	/** Obtain the ids of the automata in the network. */
	virtual automaton_id_set get_automata() const;

	/** \brief Obtain the locations satisfying the constraints lcons. */
	virtual location_id_set
			get_locations(const location_constraint_set& lcons) const;

	/** \brief Obtain the single location satisfying the complete constraints lcons. */
//	virtual location_id
//			get_location_id(const location_constraint_set& lcons) const;

	/** \brief Canonicalize the constraint con and add the resulting
	 * constraints to lcons; return true if any
	 * changes were made. */
	virtual bool
			canonicalize_location_constraint(const automaton_id& aut_id,
					const location_constraint& con,
					location_constraint_set& lcons) const;

	/* \} */
	// --------------------------------------------
	/** \name Hybrid automaton interface
	 *  \{ */
	// --------------------------------------------

	/** Creates a new hybrid automaton of the same type as h1. */
	virtual hybrid_automaton* create() const;

	/** Creates a new pairwise_hybrid_automaton_network. */
	virtual pairwise_hybrid_automaton_network* create_network() const;

	/** Add variable vid.
	 *
	 * The variable is registered as input (uncontrolled) variable if
	 * is_input is true, otherwise it is considered as controlled.
	 * It is registered as having const dynamics if is_const is true. */
	virtual void add_variable(const variable_id& vid, bool is_input = false,
			bool is_const = false);
	/** Add the variables vars, of which inp_vars are inputs (uncontrolled) and
	 * const_vars have const-dynamics.
	 *
	 * inp_vars and const_vars must be a subset of vars.
	 *
	 * The controlled variables are those in vars that are not in inp_vars.
	 */
	virtual void add_variables(const variable_id_set& vars,
			const variable_id_set& inp_vars, const variable_id_set& const_vars);
	/** Return the variables of *this (all, including input variables). */
	virtual const variable_id_set& get_variable_ids() const;
	/** Return the input variables of *this. */
	virtual const variable_id_set& get_input_variables() const;
	/** Return the const-dynamics variables of *this. */
	virtual const variable_id_set& get_const_variables() const;

	virtual const symbolic_state_collection_ptr& get_initial_states() const;

	virtual void set_initial_states(
			const symbolic_state_collection_ptr& sstate_set);

	/** Returns a pair of iterators <ibeg,iend>. ibeg points to the first outgoing transition of location l with label a,
	 * iend points just beyond the last one. The iterators can be used in a standard STL-like fashion to loop
	 * through the outgoing transitons:
	 * \code for (transition_const_iterator i=ibeg;i!=iend;++i) ...
	 * */
	virtual std::pair<transition_const_iterator, transition_const_iterator>
	get_outgoing_transitions(location_id l, label_id a) const;

	/** Compute the post-image of the transition with id trans
	 * on the locations in the discrete set *dset. */
	//virtual discrete::discrete_set::ptr post(const transition_id& trans,
	//		const discrete::discrete_set::const_ptr& dset) = 0;

	virtual transition_ptr get_transition(const transition_id& id) const;

	virtual location_ptr get_location(const location_id& id) const;
	virtual location_id get_location_id(std::string loc_name) const;
	/** Returns a pair of iterators <ibeg,iend>. ibeg points to the first location,
	 * iend points just beyond the last one. The iterators can be used in a standard STL-like fashion to loop
	 * through the locations:
	 * \code for (location_const_iterator i=ibeg;i!=iend;++i) ... */
	virtual std::pair<location_const_iterator, location_const_iterator>
			get_locations() const;

	/** Add transition and maintain the cache of outgoing transitions.
	 * The label of the transition must afterwards be in get_labels(). */
	virtual transition_id add_transition(const transition_ptr& trans, bool check_emptiness = true);

	/** Add a location.
	 *  \attention The locations should be added before their ingoing or outgoing transitions are added. */
	virtual location_id add_location(const location_ptr& loc);

	virtual const label_id_set& get_labels() const;
	/** Add the label with id lab to the alphabet of *this. */
	virtual void add_label(const label_id& lab);

	virtual void accept(hybrid_automaton_visitor& v);

	/**
	 * Output as a stream of characters.
	 */
	virtual void print(std::ostream& os) const;

	/** The set_name/get_name methods are redirected in order to
	 * have consistent naming.
	 */
	virtual const std::string& get_name() const;
	virtual void set_name(std::string s);

	/* \} */

protected:
	typedef std::map<automaton_id, hybrid_automaton::ptr>
			automaton_id_to_ptr_map;

	/** For every location of *aut_it that satisfies lcons, add the positive
	 * location constraint to fullcons and call ++aut_it.
	 * When endit is reached, obtain the location_id corresponding to fullcons and
	 * insert it into locset.
	 *
	 * Assume that lcons.is_empty() is false.
	 */
	void add_location_constraint(
			automaton_id_to_ptr_map::const_iterator aut_it,
			const automaton_id_to_ptr_map::const_iterator& endit,
			const location_constraint_set& lcons,
			location_constraint_set fullcons, location_id_set& locset) const;

	/** Set aut as the current implementation */
	void set_impl(const hybrid_automaton::ptr& aut);

	bool has_no_automata() const;

private:
	void add_base_automata(automaton_id aut_id);

	hybrid_automaton::ptr my_impl;
	automaton_id_set my_automata;
	//pw_network_node* my_root;
	automaton_id_to_ptr_map my_base_automata; // automata that are not networks
	automaton_id_to_ptr_map my_impl_automata;
};

}

#endif /* PAIRWISE_HYBRID_AUTOMATON_NETWORK_H_ */
