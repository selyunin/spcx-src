#ifndef DISCRETE_SET_H_
#define DISCRETE_SET_H_

#include <iostream>
#include "boost/shared_ptr.hpp"
#include "utility/printable.h"
#include "utility/simple_iterators/collection_base.h"
#include "utility/dispatching/double_dispatch.h"

#include "core/discrete/discrete_state.h"
#include "core/discrete/discrete_set_declarations.h"
#include "core/hybrid_automata/automaton_id.h"

namespace discrete {

/** An abstract class for representing sets of discrete states (locations).
 *
 * @todo Should the objects be pairs of (automaton_id,location_id)? I don't think so.
 */

class discrete_set: public printable,
		public simple_iterators::collection_const_base<discrete_state> {
public:
	typedef boost::shared_ptr<discrete_set> ptr;
	typedef boost::shared_ptr<const discrete_set> const_ptr;
	typedef dispatching::dispatcher<discrete_set_typelist> const_visitor;
	typedef hybrid_automata::automaton_id automaton_id;

	/** Virtual destructor */
	virtual ~discrete_set();

	/** Creates a new object of the same type as *this representing the emtpy
	 * set.
	 *
	 * Serves as a virtual constructor. */
	virtual discrete_set* create_empty() const = 0;

	/** Creates an identical copy of *this. */
	virtual discrete_set* clone() const = 0;

	/** Returns true if the set is empty and false otherwise. */
	virtual bool is_empty() const;

	/** Returns true if and only if this set is disjoint with *ds.
	 *
	 * The default implementation enumerates the objects and checks
	 * for pairwise emptyness.
	 */
	virtual bool is_disjoint_from(const object_type& loc) const;
	virtual bool is_disjoint_from(const discrete_set::const_ptr& ds) const;

	/** Returns true if and only if *this contains loc.
	 */
	virtual bool contains(const object_type& loc) const = 0;

	/** Returns true if and only if *ds is contained in *this.
	 *
	 * The default implementation enumerates the objects of ds and checks
	 * for containment (they all have to be contained).
	 */
	virtual bool contains(const discrete_set::const_ptr& ds) const;

	/** Output as a stream of characters.
	 */
	virtual void print(std::ostream& os) const = 0;

	/** Existentially quantify over aut_id
	 *
	 * The result contains no restrictions on the locations of aut_id.
	 * @note It might be necessary to canonicalize first, otherwise
	 * aut_id might not occur in the constraints.
	 */
	virtual void existentially_quantify(automaton_id aut_id) = 0;

	// --------------------------------------------

	/**  Accept a const_visitor. A const_visitor must provide the function
	 * <code> void dispatch(const T* c) </code>
	 * for all derived classes of discrete_set that are listed in
	 * discrete_set_typelist (see discrete_set_declarations). */
	virtual void accept(const_visitor& d) const = 0;
};

}

#endif /*DISCRETE_SET_H_*/
