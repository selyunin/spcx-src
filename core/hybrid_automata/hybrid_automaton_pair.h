/*
 * hybrid_automaton_pair.h
 *
 *  Created on: Aug 26, 2009
 *      Author: frehse
 */

#ifndef HYBRID_AUTOMATON_PAIR_H_
#define HYBRID_AUTOMATON_PAIR_H_

#include <boost/bimap.hpp>
#include "core/hybrid_automata/hybrid_automaton_wrapper.h"

namespace hybrid_automata {

/** Hybrid automaton implementation to represent a composition of two
 * hybrid automata. It explicitly instantiates the composition as a new
 * hybrid automaton.
 *
 * To create the composition, instantiate with init(hcomp,h1,h2).
 * The automaton hcomp is the automaton used to store locations
 * and transitions created on-the-fly by *this. During the init
 * call, *this swaps name and id with hcomp, so that from
 * the outside it now looks like *this is hcomp.
 * This is done so that composition (or other wrappers) can be stacked
 * arbitrarily, with the outside layer keeping the same id and
 * name.
 */
class hybrid_automaton_pair: public hybrid_automaton_wrapper {
public:
	typedef boost::shared_ptr<hybrid_automaton_pair> ptr;
	typedef boost::shared_ptr<const hybrid_automaton_pair> const_ptr;

	hybrid_automaton_pair();

	/** Create a hybrid automaton consisting of the parallel composition of
	 * two hybrid automata. The implementation is constructed using hcomp,
	 * which is assumed to be empty.
	 *
	 * @note The name and id of the composed automaton are those of hcomp.
	 * Internally, the wrapper automaton doing the on-the-fly composition
	 * (*this) swaps name and id with hcomp.
	 * */
	void init(hybrid_automaton::ptr hcomp, hybrid_automaton::ptr h1,
			hybrid_automaton::ptr h2);

	virtual ~hybrid_automaton_pair();

	virtual const symbolic_state_collection_ptr& get_initial_states() const;
	virtual void set_initial_states(const symbolic_state_collection_ptr& sstate_set);

	/** Returns a pair of iterators <ibeg,iend>. ibeg points to the first outgoing transition of location l with label a,
	 * iend points just beyond the last one. The iterators can be used in a standard STL-like fashion to loop
	 * through the outgoing transitons:
	 * \code for (transition_const_iterator i=ibeg;i!=iend;++i) ...
	 * */
	virtual std::pair<transition_const_iterator, transition_const_iterator>
			get_outgoing_transitions(location_id l, label_id a) const;

	/** Returns a pair of iterators <ibeg,iend>. ibeg points to the first location,
	 * iend points just beyond the last one. The iterators can be used in a standard STL-like fashion to loop
	 * through the locations:
	 * \code for (location_const_iterator i=ibeg;i!=iend;++i) ... */
	virtual std::pair<location_const_iterator, location_const_iterator> get_locations() const;

	/** \brief Obtain the locations satisfying the constraints lcons,
	 * assuming lcons constrains only the two children (no other
	 * automaton_ids).
	 *
	 * Note: there is no emptyness check in order to avoid excessive emptiness
	 * checking on recursive calls.
	 * */
	virtual location_id_set get_locations(const location_constraint_set& lcons) const;

	/** \brief Canonicalize the constraint con and add the resulting
	 * constraints to lcons; return true if any
	 * changes were made. */
	virtual bool canonicalize_location_constraint(const automaton_id& aut_id, const location_constraint& con,
			location_constraint_set& lcons) const;

	/** Checks whether the location (id1,id2) is already in the cache.
	 * If yes, store its id in loc and return true. If no, return false, leaving
	 * loc unchanged.
	 */
	bool find_location_id(location_id& loc, const location_id& id1, const location_id& id2) const;

	/** Returns the id of location (id1,id2). If the location does not yet
	 * exist, it is created (including composing time constraints etc.).
	 *
	 * @note This operation only makes sense if both left_child and right_child
	 * are defined.
	 */
	location_id get_or_add_location_id(const location_id& id1, const location_id& id2);

	/** Accept a visitor.
	 *
	 * The visitor is handed to the children first, then the call is forwarded to
	 * the wrapper, which hands it to the implementation.
	 * Calling the children first should ensure that any update on the
	 * implementation (re-composing of locations etc.) uses the visited
	 * children.
	 */
	virtual void accept(hybrid_automaton_visitor& v);


protected:
	typedef std::pair<location_id, location_id> loc_id_pair;

	/** Add the composition of the two symbolic states to res.
	 *
	 * @note res is passed as a parameter so that its implementation type is
	 * determined by the caller.
	 */
	void add_composed_states(symbolic_state_collection_ptr& res,
			const symbolic_state_const_ptr& s1, const symbolic_state_const_ptr& s2);

	/** Compute the parallel composition of the states in *s1 and *s2. */
	symbolic_state_collection_ptr compose(const symbolic_state_collection_const_ptr& s1,
			const symbolic_state_collection_const_ptr& s2);

	/** Instantiate the outgoing transitions of location l. */
	void update_outgoing_transitions(const location_id& l) const;

	/** Return the location ids of the children that correspond to location
	 * l of *this. Throws if no such location exists.
	 */
	const loc_id_pair& get_child_locations(const location_id& l) const;

private:
	hybrid_automaton::ptr left_child;
	hybrid_automaton::ptr right_child;

	typedef boost::bimap<loc_id_pair, location_id> pair_to_loc_bimap;
	pair_to_loc_bimap my_pair_to_loc_bimap;

	typedef std::map<location_id, bool> location_bool_map;
	location_bool_map my_loc_outgoing_uptodate_map;

	/** Disallow default constructor. */

};

}

#endif /* HYBRID_AUTOMATON_PAIR_H_ */
