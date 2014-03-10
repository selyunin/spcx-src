#ifndef TRANSITION_H_
#define TRANSITION_H_

#include <set>
#include <map>
//#include "../../utility/shared_ptr_user.h"
#include "core/continuous/continuous_set.h"
#include "core/discrete/discrete_set.h"
#include "core/hybrid_automata/named_label.h"
#include "core/hybrid_automata/location.h"


/** Forward declaration of classes used in header file. */
namespace continuous {
class continuous_set_transform;
typedef boost::shared_ptr<const continuous_set_transform> continuous_set_transform_const_ptr;
}

namespace hybrid_automata {

typedef std::set<label_id> label_id_set;

/** The jump constraints define, together with the invariant of the
 * target location, the continuous semantics of a transition:
 * -# a state must satisfy the guard : intersect with the guard set
 * -# the state is transformed into a set of states by the transform :
 * apply the transform
 * -# the resulting states must satisfy the invariant of the target
 * transition : intersect with target invariant
 *
 * For the sake of an easier creation process, it is allowed to define
 * jump constraints with null guard.
 * The semantics of the null guard is that it is always satisfied.
 *
 * @note Future extensions of the jump constraints might include flags such
 * as ASAP.
 */
class jump_constraints {
public:
	/** Construct without a guard. */
	explicit jump_constraints(continuous::continuous_set_transform_const_ptr t) :
		my_guard(continuous::continuous_set::const_ptr()), my_transform(t) {
	}
	;
	/** Construct with guard g, transform t. */
	jump_constraints(continuous::continuous_set::const_ptr g,
			continuous::continuous_set_transform_const_ptr t) :
		my_guard(g), my_transform(t) {
	}
	;
	const continuous::continuous_set::const_ptr& get_guard() const {
		return my_guard;
	}
	;
	void set_guard(continuous::continuous_set::const_ptr g) {
		my_guard = g;
	}
	;
	const continuous::continuous_set_transform_const_ptr& get_transform() const {
		return my_transform;
	}
	;
	void set_transform(continuous::continuous_set_transform_const_ptr t) {
		my_transform = t;
	}
	;
protected:
	/** The guard of the transition. May be a null pointer to signify that the
	 * guard is always satisfied. */
	continuous::continuous_set::const_ptr my_guard;
	/** The transform of the transition. */
	continuous::continuous_set_transform_const_ptr my_transform;
};

/** A transition class for explicit automata. It has a source and target location,
 * and a transformation of the continuous set associated to it.
 * It is a purely syntactic representation of the transition:
 * - the automaton is supposed to compute deal with transforming discrete sets
 * (as opposed to a single location).
 * - the jump_constraints describe whatever happens to the continuous set.
 * */
class transition {
public:
	typedef transition my_type;
	typedef boost::shared_ptr<my_type> ptr;
	typedef boost::shared_ptr<const my_type> const_ptr;

	virtual ~transition() {
	}
	;
	/** \brief Virtual constructor. */
	virtual transition* create(location_id source_loc, label_id l, jump_constraints jcons,
			location_id target_loc) const = 0;

	virtual label_id get_label() const = 0;
	virtual const location_id& get_source() const = 0;
	virtual const location_id& get_target() const = 0;
	virtual const jump_constraints& get_jump_constraints() const = 0;
	/** @note Modifying labels, source or target is forbidden since this would
	 * require an update of outgoing transition information in the automaton.
	 */
	/*
	virtual void set_label(label_id a) = 0;
	virtual void set_source(const location_id& l) = 0;
	virtual void set_target(const location_id& l) = 0;
	*/
	virtual void set_jump_constraints(const jump_constraints& t) = 0;
};

}

#endif /*TRANSITION_H_*/
