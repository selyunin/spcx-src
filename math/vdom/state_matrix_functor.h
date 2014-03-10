/*
 * state_matrix_functor.h
 *
 *  Created on: Jul 28, 2011
 *      Author: reymann
 */

#ifndef STATE_MATRIX_FUNCTOR_H_
#define STATE_MATRIX_FUNCTOR_H_

#include "vdom_vector.h"
#include "math/matrix.h"
#include "affine_map.h"

namespace math {

/** A class for representing function objects that map a state to a matrix of scalars.
 *
 * A state is a point in a space of Reals, where each dimension is
 * associated with a variables.*/
template<typename scalar_type>
class state_matrix_functor {
public:
	/** Member types */
	typedef math::vdom_matrix<scalar_type> matrix;

	/** A class for representing state in a space of Reals, where each dimension is
	 * associated with a variables.
	 */
	typedef vdom_vector<scalar_type> state;

	/** Virtual destructor for possible derived classes */
	virtual ~state_matrix_functor() {
	}
	;

	/** Pure virtual interface for mapping a state to another state */
	virtual matrix map(const state& x) const = 0;
};

/** A type conversion class for state vector functors
 *
 * Allows to use a functor that is defined for a different scalar type,
 * here called impl_type. The input is converted from scalar_type to impl_type,
 * and the output is converted back from impl_type to scalar_type. */
template<typename scalar_type,typename impl_type>
class convert_state_matrix_functor: public state_matrix_functor<scalar_type> {
public:
	/** Member types */
	typedef typename state_matrix_functor<scalar_type>::state state;
	typedef typename state_matrix_functor<scalar_type>::matrix matrix;
	typedef state_matrix_functor<scalar_type> impl_functor;
	typedef typename impl_functor::state impl_state;

	/** Construct a conversion functor from a given implementation functor f */
	convert_state_matrix_functor(const impl_functor& f) :
			my_f(f) {
	}
	;

	/** Virtual destructor */
	virtual ~convert_state_matrix_functor() {
	}
	;

	/** Map the state x using the implementation functor */
	virtual matrix map(const state& x) const {
		impl_state xo = x.template convert_to<impl_type> ();
		return my_f.map(xo).template convert_to<scalar_type> ();
	}
	;
private:
	const impl_functor& my_f;
};

}

#endif /* STATE_MATRIX_FUNCTOR_H_ */
