#ifndef DISCRETE_SET_IMPL_H_
#define DISCRETE_SET_IMPL_H_

#include <list>
#include <iostream>
#include "core/discrete/discrete_set.h"
#include <boost/enable_shared_from_this.hpp>

namespace discrete {

class discrete_set_stl_set;
typedef boost::shared_ptr<discrete_set_stl_set> discrete_set_stl_set_ptr;
typedef boost::shared_ptr<const discrete_set_stl_set> discrete_set_stl_set_const_ptr;

class discrete_set_stl_set : public discrete_set, public boost::enable_shared_from_this<discrete_set_stl_set>
{
public:
	typedef boost::shared_ptr<discrete_set_stl_set> ptr;
	typedef boost::shared_ptr<const discrete_set_stl_set> const_ptr;

	typedef discrete_set::object_type object_type;
	typedef std::list<object_type> container_type;
	typedef discrete_set::const_iterator const_iterator;
	virtual const const_iterator& begin() const {return get_const_begin(_myset.begin()); };
	virtual const const_iterator& end() const {return get_const_end(_myset.end()); };

	/**
	 * \brief Constructor
	 *
	 * Sets the implemention type of the discrete set ot be set.
	 */
	discrete_set_stl_set();

	/**
	 * \brief Virtual destructor
	 */
	~discrete_set_stl_set();

	/** Virtual constructor. */
	discrete_set_stl_set* create_empty() const;

	/** Deep copy. */
	discrete_set_stl_set* clone() const;

	/** Return a shared_ptr to *this.
	 */
	discrete_set_stl_set::ptr get_ptr() {
		return boost::enable_shared_from_this<discrete_set_stl_set>::shared_from_this();
	}
	;

	/** Return a shared_ptr to const *this.
	 */
	discrete_set_stl_set::const_ptr get_const_ptr() const {
		return boost::enable_shared_from_this<discrete_set_stl_set>::shared_from_this();
	}
	;

	/**
	 * \brief C++ class accessor function.
	 *
	 * \returns the set<T> member of the class.
	 */
	const container_type& get_stl_set() const;

	/**
	 * \brief Add a new element.
	 */
	void add(const object_type& e);

	/**
	 * \brief Utility function that returns the set size.
	 *
	 * \return Returns the number of elements in the set.
	 */
	std::size_t get_size() const;

	/**
	 * \brief Basic Boolean set manipulation function.
	 *
	 * Checks for the emptyness of the set.
	 * \return true if and only if the set is empty.
	 */
	bool is_empty() const;

	/** Returns true if and only if *this contains loc.
	 */
	bool contains(const object_type& loc) const;

//	/** Returns true if *dsp is disjoint with *this.
//	 */
//	bool is_disjoint_from(const discrete_set::const_ptr& dsp) const;

	/** \brief Union
	 *
	 * The union of this set and the set pointed by ds as assigned to this set.
	 * \param dsp A pointer to a discrete set.
	 */
	void union_assign(const discrete_set::const_ptr& dsp);

	/** Add loc to *this */
	void union_assign(const object_type& loc);

	/** Intersect every object of *this with loc. */
	void intersection_assign(const object_type& loc);

	/** \brief Intersection
	 *
	 * The intersection of this set and the set pointed by ds is assigned to this.
	 * \param dsp A pointer to a discrete set.
	 */
	void intersection_assign(const discrete_set_stl_set& dsp);

	/** Remove loc from *this */
	void difference_assign(const object_type& loc);

	/**
	 *  \brief Basic boolean set manipulation function.
	 *
	 * The difference of this set and the set pointed by ds is assigned to this.
	 * \param dsp A pointer to a discrete set.
	 */
	void difference_assign(const discrete_set::const_ptr& dsp);

	/**
	 * \brief Clears the set.
	 */
	void create_empty();

	/** Existentially quantify over aut_id
	 *
	 * The result contains no restrictions on the locations of aut_id.
	 * @note It might be necessary to canonicalize first, otherwise
	 * aut_id might not occur in the constraints.
	 */
	virtual void existentially_quantify(automaton_id aut_id);

	/**
	 * Output as a stream of characters.
	 */
	void print(std::ostream& os) const;

	/**  Accept a const_visitor. A const_visitor must provide the function
	 * <code> void dispatch(const T* c) </code>
	 * for all derived classes of discrete_set that are listed in
	 * discrete_set_typelist (see discrete_set_declarations). */
	void accept(const_visitor& d) const;

protected:
	/** Remove empty objects */
	void remove_empty();

private:
	container_type _myset;   // assuming locations ids will be integers.
};

}

#endif /*DISCRETE_SET_IMPL_H_*/
