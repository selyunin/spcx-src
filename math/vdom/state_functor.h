/*
 * state_functor.h
 *
 *  Created on: Jun 15, 2011
 *      Author: frehse
 */

#ifndef STATE_FUNCTOR_H_
#define STATE_FUNCTOR_H_

#include "vdom_vector.h"
#include "affine_map.h"

namespace math {

/** A class for representing function objects that map a state to another state.
 *
 * A state is a point in a space of Reals, where each dimension is
 * associated with a variables.*/
template<typename scalar_type>
class state_functor {
public:
	/** A class for representing state in a space of Reals, where each dimension is
	 * associated with a variables.
	 */
	typedef vdom_vector<scalar_type> state;

	typedef boost::shared_ptr<state_functor> ptr;
	typedef boost::shared_ptr<const state_functor> const_ptr;

	/** Virtual destructor for possible derived classes */
	virtual ~state_functor() {
	}
	;

	/** Pure virtual interface for mapping a state to another state */
	virtual state map(const state& x) const = 0;
};

}

#endif /* STATE_FUNCTOR_H_ */
