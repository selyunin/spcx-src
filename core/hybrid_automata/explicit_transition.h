/*
 * explicit_transition.h
 *
 *  Created on: Jun 22, 2009
 *      Author: frehse
 */

#ifndef EXPLICIT_TRANSITION_H_
#define EXPLICIT_TRANSITION_H_

#include "core/hybrid_automata/transition.h"

namespace hybrid_automata {

/** A transition class for explicit automata. It has a source and target location,
 * and a transformation of the continuous set associated to it.
 * It is a purely syntactic representation of the transition:
 * - the automaton is supposed to compute deal with transforming discrete sets
 * (as opposed to a single location).
 * - the jump_constraints describe whatever happens to the continuous set.
 * */
class explicit_transition: public transition {
public:
	typedef explicit_transition my_type;
	typedef boost::shared_ptr<my_type> ptr;
	typedef boost::shared_ptr<const my_type> const_ptr;

	explicit_transition(location_id source_loc, label_id l, jump_constraints jcons,
			location_id target_loc);
	explicit_transition(location_id source_loc, std::string label_name, jump_constraints jcons,
			location_id target_loc);
	virtual ~explicit_transition();

	/** \brief Virtual constructor. */
	virtual explicit_transition* create(location_id source_loc, label_id l, jump_constraints jcons,
			location_id target_loc) const;

	virtual label_id get_label() const;
	virtual const location_id& get_source() const;
	virtual const location_id& get_target() const;
	virtual const jump_constraints& get_jump_constraints() const;
	/*
	virtual void set_label(label_id a);
	virtual void set_source(const location_id& l);
	virtual void set_target(const location_id& l);
	*/
	virtual void set_jump_constraints(const jump_constraints& t);

private:
	location_id my_source_loc;
	location_id my_target_loc;
	label_id my_label;
	jump_constraints my_jump_constraints;
};

}

#endif /* EXPLICIT_TRANSITION_H_ */
