#ifndef LOCATION_CONSTRAINT_SET_H_
#define LOCATION_CONSTRAINT_SET_H_

#include <map>
#include "core/hybrid_automata/location_id.h"
#include "core/hybrid_automata/automaton_id.h"
#include "utility/printable.h"

#include "location_constraint.h"

namespace hybrid_automata {

/** A location_constraint_set represents sets of locations of a network of
 * automata.
 * For a given set of automata, the set of location_constraints represents
 * the Cartesian product of the locations that satisfy
 * the constraints.
 *
 * A location_constraint_set is called complete for a set of automata if
 * it attributes exactly one location to each of the automata.
 *
 * @note To simplify the implementation, the empty set is represented by
 * a location_constraint over aut=0.
 *
 * There can be multiple inequality constraints for a given automaton, but
 * only one equality constraint. If there is an equality constraint, then
 * there are no inequality constraints (they are implied by it). */
class location_constraint_set: public printable {
public:
	typedef std::multimap<automaton_id, location_constraint> container_impl_type;
	typedef container_impl_type::iterator iterator;
	typedef container_impl_type::const_iterator const_iterator;
	typedef std::pair<container_impl_type::iterator,container_impl_type::iterator> iterator_pair;
	typedef std::pair<container_impl_type::const_iterator,container_impl_type::const_iterator> const_iterator_pair;
	typedef container_impl_type::size_type size_type;

	/** Constructs a universe location constraint set. */
	location_constraint_set();
	explicit location_constraint_set(const automaton_id& aut_id, const location_id& id);
	virtual ~location_constraint_set();

	/** Returns the number of constraints in *this. If *this is empty, the number is 1. */
	size_type size() const;

	/** Return the constraint assigned to aut_id. If there is none,
	 * return the null pointer.
	 */
	virtual const_iterator_pair get_constraints(const automaton_id& aut_id) const;

	virtual iterator begin();
	virtual iterator end();
	virtual const_iterator begin() const;
	virtual const_iterator end() const;
	virtual bool is_empty() const;
	virtual bool is_universe() const;

	/** Returns true if the locations satisfying *this are satisfying
	 * lcons, i.e., if lcons is contained in *this.
	 *
	 * This is the case if and only if lcons is empty, or if for every
	 * constraint in *this there is a noncontradicting constraint in lcons,
	 * i.e.,
	 * 1) the constraint must be equal, or
	 * 2) the constraint in *this is inequality and in ds it's equality
	 *    for a different location.
	 */
	virtual bool contains(const location_constraint_set& lcons) const;

	/** Returns true if none of the locations satisfying *this are satisfying
	 * lcons, i.e., if their intersection is empty.
	 *
	 * This is the case if and only if either one is empty, or if they
	 * contain a contradicting constraint.
	 */
	virtual bool is_disjoint_from(const location_constraint_set& lcons) const;

	virtual void set_empty();

	/** Add the constraint aut_id=loc_id. If this contradicts an existing
	 * constraint, the result is empty. */
	virtual void add_constraint(automaton_id aut_id, location_constraint con);

	/** Add the constraints in lcons. */
	virtual void intersection_assign(const location_constraint_set& lcons);



	/** Returns the automata restricted by the constraints in *this.
	 * The set is empty if *this is empty or universal.
	 */
	virtual automaton_id_set get_automata() const;

	/** Returns true if *this defines exactly one location to
	 * every automaton in aset.
	 */
	virtual bool is_complete(const automaton_id_set& aset) const;

	/** Replaces constraints on automaton aut1 with aut2
	 *
	 * Returns true if any changes were made. */
	virtual bool map(automaton_id aut1, automaton_id aut2);

	/**
	 * Existentially quantify all elements with the first key member as aut_id.
	 *
	 * @param aut_id The automaton to be projected out of the network of automata
	 */
	virtual void existentially_quantify(automaton_id aut_id);

	/** Comparison operator */
	virtual bool operator==(const location_constraint_set& lcs) const;
	virtual bool operator!=(const location_constraint_set& lcs) const;
	/** Lexicographical comparison between constraints */
	virtual bool operator<(const location_constraint_set& lcs) const;

	/** Output as a stream of characters. */
	virtual void print(std::ostream& os) const;

protected:
	/** Return true if there is a constraint in *this that has a
	 * contradicting constraint in lcons.
	 */
	bool strictly_contradicts(const location_constraint_set& lcons) const;

	/** Return the constraint assigned to aut_id. If there is none,
	 * return the null pointer.
	 */
	virtual iterator_pair get_constraints(const automaton_id& aut_id);

private:
	container_impl_type my_container;
};

}

#endif /*LOCATION_CONSTRAINT_SET_H_*/
