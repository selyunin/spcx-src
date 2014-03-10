/*
 * hybrid_automaton_wrapper.h
 *
 *  Created on: Sep 10, 2009
 *      Author: frehse
 */

#ifndef HYBRID_AUTOMATON_WRAPPER_H_
#define HYBRID_AUTOMATON_WRAPPER_H_

#include "core/hybrid_automata/hybrid_automaton.h"

/** Forward declarations */
namespace hybrid_automata {
class symbolic_state;
typedef boost::shared_ptr<symbolic_state> symbolic_state_ptr;
typedef boost::shared_ptr<const symbolic_state> symbolic_state_const_ptr;
class location_constraint_set;
}

namespace hybrid_automata {

/** A wrapper class for hybrid automata
 *
 * This class forwards (almost all) calls to a
 * hybrid automaton implementation. Derived classes can
 * intercept these calls and add their own functionality,
 * calling upon the hybrid automaton implementation as needed.
 *
 *
 * The wrapper can be used for various purposes:
 * - adding functionality by executing code before or after passing it to
 *   the implementation (decorator pattern),
 * - changing/implementing functionality (strategy pattern),
 * - on-demand implementation of expensive operations (virtual proxy pattern),
 * - composition (composite pattern).
 *
 * @note The handling of the id is crucial: all calls to the original
 * automaton must be redirected to the wrapper, which means that
 * the wrapper must adopt the id and name of the implementation.
 * The implementation is called "...IMPLXXX", where XXX is the original id of
 * the wrapper.
 *
 * @note The set_name/get_name methods are not redirected in order to
 * have consistent naming. Since identities are swapped in the constructor,
 * this poses no problems.
 */
class hybrid_automaton_wrapper: public hybrid_automaton {
public:
	typedef boost::shared_ptr<hybrid_automaton_wrapper> ptr;
	typedef boost::shared_ptr<const hybrid_automaton_wrapper> const_ptr;

	/** Return a shared_ptr to *this. */
	ptr get_ptr();

	/** Return a shared_ptr to const *this. */
	const_ptr get_const_ptr() const;

	hybrid_automaton_wrapper();

	/** Create a wrapper for the implementation himpl. The wrapper
	 * swaps ids with himpl in order to redirect all calls. */
	void set_impl(hybrid_automaton::ptr himpl);

	virtual ~hybrid_automaton_wrapper();

	/** Tries to creates a new hybrid automaton of the same type as *this,
	 * i.e., the same implementation in the same wrapper.
	 * If this is not possible, the create is passed on to the implementation
	 * automaton of the wrapper.
	 *
	 * Derived classes should construct an object of the
	 * same class and set as implementation the create() of the
	 * implementation. That way nested wrappers are properly
	 * instantiated. */
	virtual hybrid_automaton* create() const;

	/** Add variable vid.
	 *
	 * The variable is registered as input (uncontrolled) variable if
	 * is_input is true, otherwise it is considered as controlled.
	 * It is registered as having const dynamics if is_const is true. */
	virtual void add_variable(const variable_id& vid, bool is_input=false, bool is_const=false);
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
	/** Return the const variables of *this. */
	virtual const variable_id_set& get_const_variables() const;

	virtual const symbolic_state_collection_ptr& get_initial_states() const;

	virtual void set_initial_states(const symbolic_state_collection_ptr& sstate_set);

	/** Returns a pair of iterators <ibeg,iend>. ibeg points to the first outgoing transition of location l with label a,
	 * iend points just beyond the last one. The iterators can be used in a standard STL-like fashion to loop
	 * through the outgoing transitons:
	 * \code for (transition_const_iterator i=ibeg;i!=iend;++i) ...
	 * */
	virtual std::pair<transition_const_iterator, transition_const_iterator>
			get_outgoing_transitions(location_id l, label_id a) const;

	virtual transition_ptr get_transition(const transition_id& id) const;
	virtual location_ptr get_location(const location_id& id) const;
	virtual location_id get_location_id(std::string loc_name) const;

	/** Returns a pair of iterators <ibeg,iend>. ibeg points to the first location,
	 * iend points just beyond the last one. The iterators can be used in a standard STL-like fashion to loop
	 * through the locations:
	 * \code for (location_const_iterator i=ibeg;i!=iend;++i) ... */
	virtual std::pair<location_const_iterator, location_const_iterator> get_locations() const;

	/** \brief Obtain the locations satisfying the constraints lcons.
	 *
	 * Note: there is no emptiness check on lcons in order to avoid excessive
	 * emptiness checking on nested calls.
	 * */
	virtual location_id_set get_locations(const location_constraint_set& lcons) const;

	/** \brief Canonicalize the constraint con and add the resulting
	 * constraints to lcons; return true if any
	 * changes were made. */
	virtual bool canonicalize_location_constraint(const automaton_id& aut_id, const location_constraint& con,
			location_constraint_set& lcons) const;

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


	/** Output as a stream of characters.
	 */
	virtual void print(std::ostream& os) const;

protected:
	hybrid_automaton::ptr my_impl;

private:
	/** Disallow default constructor. */
};

}

#endif /* HYBRID_AUTOMATON_WRAPPER_H_ */
