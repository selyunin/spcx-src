#ifndef PASSED_AND_WAITING_LIST_H_
#define PASSED_AND_WAITING_LIST_H_

#include <boost/shared_ptr.hpp>
#include "utility/printable.h"

/** Forward declarations */
namespace hybrid_automata {
class hybrid_automaton;
typedef boost::shared_ptr<hybrid_automaton> hybrid_automaton_ptr;
typedef boost::shared_ptr<const hybrid_automaton> hybrid_automaton_const_ptr;
class symbolic_state;
typedef boost::shared_ptr<symbolic_state> symbolic_state_ptr;
typedef boost::shared_ptr<const symbolic_state> symbolic_state_const_ptr;
class symbolic_state_collection;
typedef boost::shared_ptr<symbolic_state_collection> symbolic_state_collection_ptr;
typedef boost::shared_ptr<const symbolic_state_collection> symbolic_state_collection_const_ptr;
}

namespace hybrid_automata {

/**
 * 	Interface for passed/waiting lists.
 *
 *  The waiting list keeps track of the symbolic states on which the post (discrete and time elapse)
 *  operators needs to be applied to explore further reachable states of the automaton.
 *  The passed list keeps track of the symbolic states that have been seen already.
 *  The waiting list keeps track of the symbolic states that are newly discovered and who
 *  need to be explored with post operators.
 *
 *	Every state on the waiting list needs to eventually end up on the passed list, possibly
 *	as a subset of another state.
 **/
class passed_and_waiting_list: public printable {
public:
	typedef boost::shared_ptr<passed_and_waiting_list> ptr;
	typedef boost::shared_ptr<const passed_and_waiting_list> const_ptr;

	// --------------------------------------------
	// Constructors
	// --------------------------------------------
	/** \name Constructors
	 *  \{ */

	/** Constructor */
	passed_and_waiting_list() {
	}
	;

	/** Virtual Constructor */
	virtual passed_and_waiting_list* create() const = 0;

	/** Virtual Deep Copy Constructor */
	virtual passed_and_waiting_list* clone() const = 0;

	/** Destructor */
	virtual ~passed_and_waiting_list() {
	}
	;

	/* \} */
	// --------------------------------------------
	/** \name Accessor and test methods
	 *  \{ */
	// --------------------------------------------

	/** Accessor method for passed list. */
	virtual symbolic_state_collection_ptr get_passed_list() const = 0;

	/** Accessor method for waiting list. */
	virtual symbolic_state_collection_ptr get_waiting_list() const = 0;

	/** Returns whether the waiting list is empty. */
	virtual bool is_waiting_list_empty() const = 0;

	/** Information about the PLWL */
	struct info {
		unsigned int pl_size;
		unsigned int wl_size;
	};

	/** Returns information about the PLWL. */
	virtual info get_info() const = 0;

	/* \} */
	// --------------------------------------------
	/** \name Basic definition methods
	 *  \{ */
	// --------------------------------------------

	/** Clear the passed list. */
	virtual void clear_passed_list() = 0;

	/** Clear the waiting list. */
	virtual void clear_waiting_list() = 0;

	/** Define a new passed list.  */
	virtual void set_passed_list(const symbolic_state_collection_ptr& new_pl) = 0;

	/** Define a new waiting list. */
	virtual void set_waiting_list(const symbolic_state_collection_ptr& new_wl) = 0;

	/* \} */
	// --------------------------------------------
	/** \name Main semantic methods
	 *  \{ */
	// --------------------------------------------

	/**
	 * Remove a state from the waiting list and hand it over for processing (adopt).
	 */
	virtual symbolic_state_ptr pick() = 0;

	/**
	 * \brief Add the symbolic state referred to by \p sstate to the passed list.
	 * \param sstate A pointer to the symbolic state to be added.
	 * \param H The hybrid automaton whose reachable set computation is desired.
	 */
	virtual void add_passed(symbolic_state_ptr sstate, const hybrid_automaton_ptr H) = 0;

	/**
	 * \brief Add the symbolic state referred to by \p sstate to the passed list.
	 * \param sstate_set A pointer to the symbolic state to be added.
	 * \param H The hybrid automaton whose reachable set computation is desired.
	 */
	virtual void add_passed(symbolic_state_collection_ptr sstate_set, const hybrid_automaton_ptr H);

	/**
	 * \brief Add the symbolic state referred to by \p sstate to the waiting list.
	 * \param sstate A pointer to the symbolic state to be added.
	 * \param H The hybrid automaton whose reachable set computation is desired.
	 */
	virtual void add_waiting(symbolic_state_ptr sstate, const hybrid_automaton_ptr H) = 0;

	/**
	 * \brief Add the symbolic state referred to by \p sstate to the waiting list.
	 * \param sstate_set A pointer to the symbolic state to be added.
	 * \param H The hybrid automaton whose reachable set computation is desired.
	 */
	virtual void add_waiting(symbolic_state_collection_ptr sstate_set, const hybrid_automaton_ptr H);
	/* \} */

	/** Output as a stream of characters.
	 */
	virtual void print(std::ostream& os) const;

	/** A generic interface to set options.
	 *
	 * The default behavior is to ignore the option. */
	virtual void set_option(std::string opt);
};

}

#endif /*PASSED_AND_WAITING_LIST_H_*/
