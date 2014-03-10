#ifndef SYMBOLIC_STATE_COLLECTION_STL_LIST_H_
#define SYMBOLIC_STATE_COLLECTION_STL_LIST_H_

#include <iostream>
#include <list>
#include "boost/shared_ptr.hpp"
#include "core/symbolic_states/symbolic_state_collection.h"

namespace hybrid_automata {

class symbolic_state_collection_stl_list;
typedef boost::shared_ptr<symbolic_state_collection_stl_list> symbolic_state_collection_stl_list_ptr;
typedef boost::shared_ptr<const symbolic_state_collection_stl_list> symbolic_state_collection_stl_list_const_ptr;


class symbolic_state_collection_stl_list :
	public symbolic_state_collection {
private:
	std::list<symbolic_state_ptr> mysymbolic_state_list;

public:
	typedef boost::shared_ptr<symbolic_state_collection_stl_list> ptr;
	typedef boost::shared_ptr<const symbolic_state_collection_stl_list> const_ptr;

	/**
	 * \brief Virtual Destructor. It's a C++ specific function.
	 */
	virtual ~symbolic_state_collection_stl_list();

	const const_iterator& begin() const {
		return get_const_begin(mysymbolic_state_list.begin());
	}
	;
	const const_iterator& end() const {
		return get_const_end(mysymbolic_state_list.end());
	}
	;


	symbolic_state_collection_stl_list();

	/** Virtual constructor. */
	symbolic_state_collection_stl_list* create() const;

	/** Virtual copy constructor. */
	symbolic_state_collection_stl_list* clone() const;

	/**
	 * \brief Utility function.
	 *
	 * \return The memory used by the object.
	 */
	std::size_t get_memory() const;

	/** \brief Makes the collection empty.
	 */
	void clear() {
		mysymbolic_state_list.clear();
	}

	/**
	 * \brief Class member Accesor function to access the STL list of symbolic states.
	 */
	const std::list<symbolic_state_ptr> get_list() const {
		return mysymbolic_state_list;
	}

	/**
	 * \brief Class member mutator function.
	 *
	 * Sets the STL list of symbolic states of \p *this to the passed parameter \p l
	 */
	void set_list(std::list<symbolic_state_ptr> l) // pass by value
	{
		mysymbolic_state_list = l; // copy assignment
	}

	/**
	 * Checks for the emptyness of *this list of symbolic state collection.
	 * \return \p true if and only if list of symbolic states is empty.
	 */
	bool is_empty() const;

	/**
	 * Returns the size of *this.
	 */
	virtual unsigned int size() const;

	/**
	 * \brief Remove the symbolic state at the front of the list and return it.
	 *
	 */
	symbolic_state_ptr erase_front();

	/**
	 * \brief Auxiliary set manipulation function.
	 *
	 * \return \p true if and only if the set pointer why the passed parameter is contained in
	 * the set pointed by this.
	 */
	bool contains(const symbolic_state_collection::ptr& sstate_set) const;

	/*! Swaps \p *this and \p s.
	 * A utility function.
	 */
	void swap(symbolic_state_collection::ptr sstate_set);
	/*! \brief Basic boolean set manipulation function.
	 *
	 * Assigns to \p *this the union of \p *this and \p *ps, and does some simplification
	 * of the list without affecting the union result.
	 */
	void union_assign(const symbolic_state_collection::ptr& sstate_set);

	/*! \brief
	 * Basic boolean set manipulation function.
	 * Assigns to \p *this the states of \p *this that are not in \p sstate
	 */
	void difference_assign(const symbolic_state_const_ptr& sstate);

	/*! \brief
	 * Basic boolean set manipulation function.
	 * Assigns to \p *this the states of \p *this that are not in \p sstate_set
	 */
	void difference_assign(const symbolic_state_collection::const_ptr& sstate_set);

	/** Adapt the discrete sets in *this according to v.
	 * Return true if adaptation successful. */
	virtual bool accept(adapt_discrete_set_visitor& v);

	/** Adapt the continuous sets in *this according to v.
	 * Return true if adaptation successful. */
	virtual bool accept(adapt_continuous_set_visitor& v);

	/*! \brief Basic boolean set manipulation function.
	 *
	 * Assigns to \p *this the union of \p *this and \p *ps
	 */
	void intersection_assign(const symbolic_state_collection::ptr& sstate_set);

	/** \brief Assigns to \p *this the states of \p *this that are not in \p sstate.
	 *
	 * This implementation is efficient than difference_assign but trades off with accuracy.
	 * The result is an overapproximation of the actual set.
	 */
	void cheap_difference_assign(const symbolic_state_ptr& sstate);

	/** \brief An efficient overapproximation of set difference.
	 *
	 * This implementation is more efficient than difference_assign but trades off with accuracy.
	 * The result is an overapproximation of the actual set.
	 */
	void cheap_difference_assign(const symbolic_state_collection::ptr& sstate_set);

	/**
	 * \brief Add a copy of \p sstate to \p *this.
	 * Does not check if \p sstate is empty.
	 */
	void add(const symbolic_state_ptr& sstate);

	/** Return an overapproximation of the states that are not already in \p *this.
	 */
	void compute_new(const symbolic_state_ptr& sstate);

	/** Add the state to \p *this and return an overapproximation
	 * of the states that were not already in \p *this.
	 */
	void add_and_return_new(const symbolic_state_ptr& sstate);

	/** Add the state to \p *this and return an overapproximation
	 * of the states that were not already in \p *this.
	 * If a state in *this is contained in sstate, it is replaced by sstate.
	 */
	void add_and_return_new_with_merging(const symbolic_state_ptr& sstate);

};

}

#endif /*SYMBOLIC_STATE_COLLECTION_IMPL_H_*/
