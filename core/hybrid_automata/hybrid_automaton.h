#ifndef HYBRID_AUTOMATON_H_
#define HYBRID_AUTOMATON_H_

#include <list>
#include <map>
#include <boost/enable_shared_from_this.hpp>
#include "math/vdom/variable.h"
#include "utility/shared_ptr_output.h"
#include "utility/printable.h"
#include "utility/simple_iterators/collection_base.h"
//#include "location.h"
//#include "transition.h"
//#include "../symbolic_states/symbolic_state_collection.h"
//#include "../discrete/singleton_set.h"

#include "core/hybrid_automata/automaton_id.h"
#include "core/hybrid_automata/named_label.h"
#include "core/hybrid_automata/location_id.h"
#include "core/hybrid_automata/transition_id.h"

/** Forward declarations */
namespace hybrid_automata {
class location_constraint;
class location_constraint_set;
class symbolic_state_collection;
typedef boost::shared_ptr<symbolic_state_collection> symbolic_state_collection_ptr;
typedef boost::shared_ptr<const symbolic_state_collection> symbolic_state_collection_const_ptr;
class transition;
typedef boost::shared_ptr<transition> transition_ptr;
typedef boost::shared_ptr<const transition> transition_const_ptr;
class jump_constraints;
class location;
typedef boost::shared_ptr<location> location_ptr;
typedef boost::shared_ptr<const location> location_const_ptr;
}
namespace discrete {
class discrete_set;
typedef boost::shared_ptr<discrete_set> discrete_set_ptr;
typedef boost::shared_ptr<const discrete_set> discrete_set_const_ptr;
}
namespace continuous {
class continuous_set;
typedef boost::shared_ptr<continuous_set> continuous_set_ptr;
typedef boost::shared_ptr<const continuous_set> continuous_set_const_ptr;
}

namespace hybrid_automata {

class hybrid_automaton_visitor;

/** Hybrid automaton class in which locations and transitions are explicitly present as objects.
 *
 * \note These objects are identified via ids, which is supposed to make them referenceable
 * independently of cloning etc. E.g., when a state set R contains location id 1 of automaton 1,
 * then this set should also be able to refer to copies or modifications of automaton 1.
 *
 * \note The label_ids are globally defined.
 *
 * \note Every transition needs a label. A transition without label (and therefore not synchronizing
 * with any other transition) is modeled with a special label called silent label.
 * Its id is named_label::silent_id(). If a transition with this label exists,
 * the silent id needs to be part of the alphabet of the automaton.
 */
class hybrid_automaton: public virtual printable, public boost::enable_shared_from_this<
		hybrid_automaton> {
public:
	typedef boost::shared_ptr<hybrid_automaton> ptr;
	typedef boost::shared_ptr<const hybrid_automaton> const_ptr;

	typedef simple_iterators::collection_const_base<transition_id>::const_iterator
			transition_const_iterator;

	typedef simple_iterators::collection_const_base<location_id>::const_iterator
			location_const_iterator;

	/** Return a shared_ptr to *this.*/
	ptr get_ptr();
	/** Return a shared_ptr to const *this.*/
	const_ptr get_const_ptr() const;

	hybrid_automaton();
	hybrid_automaton(std::string new_name);
	virtual ~hybrid_automaton();
	virtual hybrid_automaton* create() const = 0;

	/** Add variable vid.
	 *
	 * The variable is registered as input (uncontrolled) variable if
	 * is_input is true, otherwise it is considered as controlled.
	 * It is registered as having const dynamics if is_const is true. */
	virtual void add_variable(const variable_id& vid, bool is_input = false, bool is_const = false) = 0;
	/** Add the variables vars, of which inp_vars are inputs (uncontrolled) and
	 * const_vars have const-dynamics.
     *
     * inp_vars and const_vars must be a subset of vars.
	 *
	 * The controlled variables are those in vars that are not in inp_vars.
	 */
	virtual void add_variables(const variable_id_set& vars, const variable_id_set& inp_vars, const variable_id_set& const_vars);
	/** Return the variables of *this (all, including input variables). */
	virtual const variable_id_set& get_variable_ids() const = 0;
	/** Return the input variables of *this.
	 *
	 * The input variables are the variables that are not controlled. */
	virtual const variable_id_set& get_input_variables() const = 0;
	/** Return the controlled (non-input) variables of *this. */
	virtual variable_id_set get_controlled_variables() const;
	/** Return the const-dynamics variables of *this.
	 *
	 * The const-dynamics variables by definition do not change their value
	 * during any execution of the automaton.
	 * In literature, these are also called "parameters" of the automaton.
	 * */
	virtual const variable_id_set& get_const_variables() const = 0;

	virtual const symbolic_state_collection_ptr & get_initial_states() const = 0;

	virtual void set_initial_states(const symbolic_state_collection_ptr& sstate_set) = 0;

	/** Returns a pair of iterators <ibeg,iend>. ibeg points to the first outgoing transition of location l with label a,
	 * iend points just beyond the last one. The iterators can be used in a standard STL-like fashion to loop
	 * through the outgoing transitions:
	 * \code for (transition_const_iterator i=ibeg;i!=iend;++i) ...
	 * */
	virtual std::pair<transition_const_iterator, transition_const_iterator>
	get_outgoing_transitions(location_id l, label_id a) const = 0;

	/** Returns a list of pointers to outgoing transitions of the discrete set *d.
	 *
	 * The default implementation enumerates the locations of d and joins the
	 * respective iterator ranges. */
	virtual std::list<transition_id> get_outgoing_transitions(
			const discrete::discrete_set_const_ptr d, label_id a) const;

	/** Returns a list of pointers to discrete sets. For each of the
	 * discrete sets, the locations have identical time constraints
	 * with respect to the continuous set *cset,
	 * so the continuous post operator can be applied in unison.
	 * The discrete sets cover the discrete set *d.
	 * The returned sets may be a partition of *d, but doesn't have to be
	 * (they can overlap).
	 *
	 * The default implementation enumerates the locations of d and returns
	 * them as singleton sets, i.e., it considers none of them equivalent. */
	virtual std::list<discrete::discrete_set_ptr> get_time_equiv_locations(
			const discrete::discrete_set_const_ptr& d,
			const continuous::continuous_set_const_ptr& cset) const;

	/** \brief Obtain the ids of all locations satisfying the constraints lcons.
	 *
	 * The set is either empty (lcons is unsatisfiable), contains exactly one location
	 * (lcons contains a positive constraint on *this), contains more than one location
	 * (lcons contains one or more negative constraints on *this), or all locations of *this
	 * (lcons contains no constraint on *this). */
	virtual location_id_set get_locations(const location_constraint_set& lcons) const = 0;

	/** \brief Canonicalize the constraint con and add the resulting
	 * constraints to lcons; return true if any
	 * changes were made.
	 *
	 * If *this is atomic (non-network), the constraint id=con is added to lcons.
	 * Otherwise the constraint is expanded into constraints on atomic automata.
	 * The reason why id=con is added instead of get_id()=con is so that wrappers
	 * can substitute their own id. That way when using the constraint all
	 * calls are directed to the wrapper and not the implementation.
	 * Wrapped automata cannot detect the id of their wrappers, so this id
	 * has to be passed explicitly.
	 *
	 * @note A set of location constraints is canonic if it refers only to atomic
	 * (non-network) automata. An atomic automaton simply returns con in lcons.
	 * A network automaton (recursively) replaces con, if it refers
	 * to a composition, by a set of constraints on atomic automata.
	 * The set of constraints that is being introduced
	 * is a conjunction, so the result is again a location_constraint_set. */
	virtual bool canonicalize_location_constraint(const automaton_id& id, const location_constraint& con,
			location_constraint_set& lcons) const = 0;

	/** Compute the post-image of the transition with id trans
	 * on the locations in the discrete set *dset. */
	//virtual discrete::discrete_set::ptr post(const transition_id& trans,
	//		const discrete::discrete_set::const_ptr& dset) = 0;

	virtual transition_ptr get_transition(const transition_id& id) const = 0;

	/** Functions to access the transition information directly,
	 * can be overloaded by the hybrid automaton implementation with
	 * something more efficient. */
	virtual label_id get_label(const transition_id& id) const;
	virtual const location_id& get_source(const transition_id& id) const;
	virtual const location_id& get_target(const transition_id& id) const;
	virtual const jump_constraints& get_jump_constraints(const transition_id& id) const;

	virtual location_ptr get_location(const location_id& id) const = 0;
	virtual location_id get_location_id(std::string loc_name) const = 0;
	/** Returns a pair of iterators <ibeg,iend>. ibeg points to the first location,
	 * iend points just beyond the last one. The iterators can be used in a standard STL-like fashion to loop
	 * through the locations:
	 * \code for (location_const_iterator i=ibeg;i!=iend;++i) ... */
	virtual std::pair<location_const_iterator, location_const_iterator> get_locations() const = 0;

	/** Add transition and maintain the cache of outgoing transitions.
	 * The label of the transition must afterwards be in get_labels().
	 *
	 * If check_emptiness = true, emptiness checks can be performed on
	 * target invariant, guards etc. At the parse stage these objects
	 * might not implement such methods, so parser functions call
	 * add_transition with check_emptiness = false.
	 * If the transition is not added due to emptiness checks, the
	 * transition_id 0 is returned. */
	virtual transition_id add_transition(const transition_ptr& trans, bool check_emptiness = true) = 0;

	/** Add a location.
	 *  \attention The locations should be added before their ingoing or outgoing transitions are added. */
	virtual location_id add_location(const location_ptr& loc) = 0;

	virtual const label_id_set& get_labels() const = 0;
	/** Add the label with id lab to the alphabet of *this.
	 *
	 * @note The silent label should never be added. */
	virtual void add_label(const label_id& lab) = 0;

	/** Add the labels in labs to the alphabet of *this. */
	virtual void add_labels(const label_id_set& labs);

	/** Add the label with name lab_name to the alphabet of *this.
	 *
	 * @note The silent label should never be added. */
	virtual void add_label(const std::string& lab_name);

	virtual void accept(hybrid_automaton_visitor& v) = 0;

	/** Output as a stream of characters. */
	virtual void print(std::ostream& os) const = 0;

	virtual const std::string& get_name() const;
	virtual void set_name(std::string s);
	virtual const automaton_id& get_id() const;

private:
	friend class hybrid_automaton_cache; // needs direct access to id and name in order to swap them
	friend class hybrid_automaton_wrapper;
	automaton_id get_new_id();

	automaton_id my_id;
	std::string my_name;
};

}

#include "core/hybrid_automata/hybrid_automaton_visitor.h"

#endif /*HYBRID_AUTOMATON_H_*/
