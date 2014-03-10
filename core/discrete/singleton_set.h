/*
 * singleton_set.h
 *
 *  Created on: Jun 22, 2009
 *      Author: frehse
 */

#ifndef SINGLETON_SET_H_
#define SINGLETON_SET_H_

#include <boost/enable_shared_from_this.hpp>
#include "core/discrete/discrete_set.h"

namespace discrete {

/** A singleton_set is a set with exactly one element. The element can represent
 * the empty set, the universe or a single location.
 *
 * This class provides the interface between implementations based on having
 * a single location and those based on sets of locations.
 */
class singleton_set: public discrete_set,
		public boost::enable_shared_from_this<singleton_set>,
		public simple_iterators::collection_nonconst_base<discrete_set::object_type> {
public:
	typedef boost::shared_ptr<singleton_set> ptr;
	typedef boost::shared_ptr<const singleton_set> const_ptr;

	/** The kind encodes whether the set is he empty set, the universe or a
	 * single location. */
	enum kind {
		EMPTY, UNIVERSE, SINGLE
	};

	/** Constructs a universe set. */
	singleton_set();
	singleton_set(const object_type& obj);

	/** Virtual destructor */
	virtual ~singleton_set();

	/** Return a shared_ptr to *this. */
	discrete_set::ptr get_ptr();
	/** Return a shared_ptr to const *this. */
	discrete_set::const_ptr get_const_ptr() const;

	/** Virtual constructor. */
	virtual singleton_set* create_empty() const;

	/** Creates an identical copy of *this. */
	virtual singleton_set* clone() const;
	virtual const const_iterator& begin() const;
	virtual const const_iterator& end() const;

	virtual const iterator& begin();
	virtual const iterator& end();

	const object_type& get_object() const;

	/** Returns true if the set is empty and false otherwise. */
	virtual bool is_empty() const;
	/** Returns true if the set is universe and false otherwise. */
	virtual bool is_universe() const;

	/** Returns true if and only if this set is disjoint with ds. */
	virtual bool is_disjoint_from(const singleton_set& ds) const;

	/** Returns true if and only if *this contains loc. */
	virtual bool contains(const object_type& loc) const;

	/** Returns true if and only if ds is contained in *this. */
	virtual bool contains(const singleton_set& ds) const;

	/** Set *this to be the empty set. */
	virtual void set_empty();

	/** Set *this to be the universe set. */
	virtual void set_universe();

	/** Set *this to be the SINGLE set containing loc.s */
	virtual void set_single(const object_type& loc);

	/** The intersection of this set and the set pointed by ds is assigned
	 * to this.
	 */
	virtual void intersection_assign(const singleton_set& ds);

	/** Existentially quantify over aut_id
	 *
	 * The result contains no restrictions on the locations of aut_id.
	 * @note It might be necessary to canonicalize first, otherwise
	 * aut_id might not occur in the constraints.
	 */
	virtual void existentially_quantify(automaton_id aut_id);

	/** Output as a stream of characters.
	 */
	virtual void print(std::ostream& os) const;

	/**  Accept a const_visitor. A const_visitor must provide the function
	 * <code> void dispatch(const T* c) </code>
	 * for all derived classes of discrete_set that are listed in
	 * discrete_set_typelist (see discrete_set_declarations). */
	virtual void accept(const_visitor& d) const;

private:
	kind my_kind;
	object_type my_location;
};

}

#endif /* SINGLETON_SET_H_ */
