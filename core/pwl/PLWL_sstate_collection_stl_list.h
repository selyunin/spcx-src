#ifndef PLWL_sstate_collection_stl_list_H_
#define PLWL_sstate_collection_stl_list_H_

#include "core/pwl/passed_and_waiting_list.h"

/** Forward declarations */
namespace hybrid_automata {
class symbolic_state_collection_stl_list;
typedef boost::shared_ptr<symbolic_state_collection_stl_list> symbolic_state_collection_stl_list_ptr;
typedef boost::shared_ptr<const symbolic_state_collection_stl_list> symbolic_state_collection_stl_list_const_ptr;
}

namespace hybrid_automata {

class PLWL_sstate_collection_stl_list : public passed_and_waiting_list
{
private:
	symbolic_state_collection_stl_list_ptr passed_list;
										 /* List of symbolic states which are discovered to be reachable and yields no
										  * symbolic state not part of the current list on application of post
										  * operators.
										  */
	symbolic_state_collection_stl_list_ptr waiting_list;
									    /* List of symbolic states discovered to be reachable and can yield new
	 									 * reachable symbolic states not in the passed list upon application of
	 									 * post operators.
	 									 */
public:
	// --------------------------------------------
	// Constructors
	// --------------------------------------------
	/** \name Constructors
	 *  \{ */

	/** Constructor */
	PLWL_sstate_collection_stl_list();

	/** Virtual Constructor */
	virtual PLWL_sstate_collection_stl_list* create() const;

	/** Virtual Deep Copy Constructor */
	virtual PLWL_sstate_collection_stl_list* clone() const;

	/** Destructor */
	virtual ~PLWL_sstate_collection_stl_list();

	/* \} */
	// --------------------------------------------
	/** \name Accessor and test methods
	 *  \{ */
	// --------------------------------------------

	/** Accessor method for passed list. */
	symbolic_state_collection_ptr get_passed_list() const;

	/** Accessor method for waiting list. */
	symbolic_state_collection_ptr get_waiting_list() const;

	/** Returns whether the waiting list is empty. */
	bool is_waiting_list_empty() const;

	/** Returns information about the PLWL. */
	info get_info() const;

	/* \} */
	// --------------------------------------------
	/** \name Basic definition methods
	 *  \{ */
	// --------------------------------------------

	/** Clear the passed list. */
	void clear_passed_list();

	/** Clear the waiting list. */
	void clear_waiting_list();

	/** Define a new passed list.  */
	void set_passed_list(const symbolic_state_collection_ptr& new_pl);

	/** Define a new waiting list. */
	void set_waiting_list(const symbolic_state_collection_ptr& new_wl);

	/* \} */
	// --------------------------------------------
	/** \name Main semantic methods
	 *  \{ */
	// --------------------------------------------

	/**
	 * Remove a state from the waiting list and hand it over for processing (adopt).
	 */

	symbolic_state_ptr pick();

	/**
	 * \brief Add the symbolic state referred to by \p sstate to the passed list.
	 * \param sstate A pointer to the symbolic state to be added.
	 * \param H The hybrid automaton whose reachable set computation is desired.
	 */
	void add_passed(symbolic_state_ptr sstate, const hybrid_automaton_ptr H);

	/**
	 * \brief Add the symbolic state referred to by \p sstate to the waiting list.
	 * \param sstate A pointer to the symbolic state to be added.
	 * \param H The hybrid automaton whose reachable set computation is desired.
	 */
	void add_waiting(symbolic_state_ptr sstate, const hybrid_automaton_ptr H);

	/**
	 * \brief Empties the contents of the passed list.
	 *
	 */
	void clear_pl();

	/**
	 * \brief Boolean check of the emptiness of the waiting list.
	 *
     * \return <code>true</code> if and only if the waiting list is empty.
	 */
	bool is_wl_empty() const;

	/** Set an option.
	 *
	 * The following options are accepted:
	 * - MERGEP_ON (OFF) : merge passed with new states
	 */
	virtual void set_option(std::string opt);

private:
	bool my_merge_passed_with_new;
};

}

#endif /*PLWL_sstate_collection_stl_list_H_*/
