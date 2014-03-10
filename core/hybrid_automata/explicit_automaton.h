#ifndef explicit_automaton_H_
#define explicit_automaton_H_

#include "utility/id_to_shared_ptr_cache.h"
#include "core/hybrid_automata/hybrid_automaton.h"
//#include "../set_implementations/stl_set_implementations/discrete_set_stl_set.h"
//#include "explicit_transition.h"

namespace hybrid_automata {

/** Hybrid automaton implementation in which locations and transitions are stored
 * in caches using their ids.
 */
class explicit_automaton: public hybrid_automaton {
public:
	typedef id_to_shared_ptr_cache<location> location_collection;
	typedef std::map<std::string, location_id> location_name_to_id_map;
	typedef id_to_shared_ptr_cache<transition> transition_collection;
	typedef std::pair<location_id, label_id> loc_and_label_pair;
	typedef std::multimap<loc_and_label_pair, transition_id> loc_and_label_to_transition_id_map;

	typedef boost::shared_ptr<explicit_automaton> ptr;
	typedef boost::shared_ptr<const explicit_automaton> const_ptr;

	explicit_automaton();
	explicit_automaton(std::string new_name);
	virtual ~explicit_automaton();

	virtual explicit_automaton* create() const;

	virtual const symbolic_state_collection_ptr & get_initial_states() const;

	/** Returns a pair of iterators <ibeg,iend>. ibeg points to the first outgoing transition of location l with label a,
	 * iend points just beyond the last one. The iterators can be used in a standard STL-like fashion to loop
	 * through the outgoing transitons:
	 * \code for (transition_const_iterator i=ibeg;i!=iend;++i) ...
	 * */
	virtual std::pair<transition_const_iterator, transition_const_iterator>
			get_outgoing_transitions(location_id l, label_id a) const;

	/** Compute the post-image of the transition with id trans
	 * on the locations in the discrete set *dset. */
//	virtual discrete::discrete_set_ptr post(const transition_id& trans,
//			const discrete::discrete_set_const_ptr& dset);

	/** Add variable vid, and register it as input variable if is_input==true. */
	virtual void add_variable(const variable_id& vid, bool is_input=false, bool is_const = false);
	/** Return the variables of *this (all, including input variables). */
	virtual const variable_id_set& get_variable_ids() const;
	/** Return the input variables of *this. */
	virtual const variable_id_set& get_input_variables() const;
	/** Return the const-dynamics variables of *this.
	 *
	 * The const-dynamics variables by definition do not change their value
	 * during any execution of the automaton.
	 * In literature, these are also called "parameters" of the automaton.
	 * */
	virtual const variable_id_set& get_const_variables() const;

	/** Import the add_label functions from the base class **/
	using hybrid_automaton::add_label;

	/** Add the label with id lab to the alphabet of *this.
	 *
	 * @note The silent label should never be added. */
	virtual void add_label(const label_id& lab);

	/** Add a location.
	 *  \attention The locations should be added before their ingoing or outgoing transitions are added. */
	virtual location_id add_location(const location_ptr& loc);

	/** Add transition and maintain the cache of outgoing transitions.
	 *
	 * If check_emptiness = true or omitted, the transition guard and
	 * invariant of the the target location are checked for
	 * emptiness. If either is empty, the transition is not added and 0 is returned. */
	virtual transition_id add_transition(const transition_ptr& trans, bool check_emptiness = true);

	virtual const label_id_set& get_labels() const;

	virtual transition_ptr get_transition(const transition_id& id) const;
	virtual location_ptr get_location(const location_id& id) const;
	virtual location_id get_location_id(std::string loc_name) const;
	/** Returns a pair of iterators <ibeg,iend>. ibeg points to the first location,
	 * iend points just beyond the last one. The iterators can be used in a standard STL-like fashion to loop
	 * through the locations:
	 * \code for (location_const_iterator i=ibeg;i!=iend;++i) ... */
	virtual std::pair<location_const_iterator, location_const_iterator> get_locations() const;

	/** \brief Obtain the locations satisfying the constraints lcons. */
	virtual location_id_set get_locations(const location_constraint_set& lcons) const;

	/** \brief Canonicalize the constraint con and add the resulting
	 * constraints to lcons; return true if any
	 * changes were made. */
	virtual bool canonicalize_location_constraint(const automaton_id& aut_id, const location_constraint& con,
			location_constraint_set& lcons) const;

	virtual void set_initial_states(const symbolic_state_collection_ptr& sstate_set);

	virtual void accept(hybrid_automaton_visitor& v);

	/**
	 * Output as a stream of characters.
	 */
	virtual void print(std::ostream& os) const;

private:
	location_collection my_locations;
	label_id_set my_labels;
	transition_collection my_transitions;
	symbolic_state_collection_ptr my_initial_states;
	loc_and_label_to_transition_id_map my_outgoing_transitions;
	location_name_to_id_map my_location_name_to_id_map;

	variable_id_set my_variables;
	variable_id_set my_input_variables;
	variable_id_set my_const_variables;
};

}

#endif /*explicit_automaton_H_*/
