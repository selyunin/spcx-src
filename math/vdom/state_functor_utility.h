/*
 * state_functor_utility.h
 *
 *  Created on: Jun 15, 2011
 *      Author: frehse
 */

#ifndef STATE_FUNCTOR_UTILITY_H_
#define STATE_FUNCTOR_UTILITY_H_

#include "state_functor.h"

namespace math {

/** An state functor implementation for affine maps. */
template<typename scalar_type>
class affine_state_functor: public state_functor<scalar_type> {
public:
	/** Member types */
	typedef typename state_functor<scalar_type>::state state;

	/** Construct a state functor corresponding to an affine map */
	affine_state_functor(const affine_map<scalar_type>& M) :
		my_map(M) {
	}
	;

	/** Virtual destructor for possible derived classes */
	virtual ~affine_state_functor() {
	}
	;

	/** Map the state x to another state according to the function
	 * defined by *this. */
	virtual state map(const state& x) const {
		return my_map.map(x);
	}
	;
private:
	const affine_map<scalar_type>& my_map;
};

/** A type conversion class for state functors
 *
 * Allows to use a functor that is defined for a different scalar type,
 * here called impl_type. The input is converted from scalar_type to impl_type,
 * and the output is converted back from impl_type to scalar_type. */
template<typename scalar_type,typename impl_type>
class convert_state_functor: public state_functor<scalar_type> {
public:
	/** Member types */
	typedef typename state_functor<scalar_type>::state state;
	typedef state_functor<scalar_type> impl_functor;
	typedef typename impl_functor::state impl_state;

	/** Construct a conversion functor from a given implementation functor f */
	convert_state_functor(const impl_functor& f) :
		my_f(f) {
	}
	;

	/** Virtual destructor */
	virtual ~convert_state_functor() {
	}
	;

	/** Map the state x using the implementation functor */
	virtual state map(const state& x) const {
		impl_state xo = x.template convert_to<impl_type> ();
		return my_f.map(xo).template convert_to<scalar_type> ();
	}
	;
private:
	const impl_functor& my_f;
};

}

#endif /* STATE_FUNCTOR_UTILITY_H_ */
