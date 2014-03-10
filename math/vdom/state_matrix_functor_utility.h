/*
 * state_matrix_functor.h
 *
 *  Created on: Jul 28, 2011
 *      Author: reymann
 */

#ifndef STATE_MATRIX_FUNCTOR_UTILITY_H_
#define STATE_MATRIX_FUNCTOR_UTILITY_H_

#include "state_matrix_functor.h"
#include "state_functor.h"
#include "math/ode_solving/traj_simu/lin_constraint_evaluation.h"

namespace math {

/** A functor that produces the same matrix, whatever inputs it takes
 * (used for Jacobian of affine dynamics)
 */
template<typename scalar_type>
class constant_state_matrix_functor: public state_matrix_functor<scalar_type> {
public:
	/** Member types */
	typedef typename state_matrix_functor<scalar_type>::state state;
	typedef typename state_matrix_functor<scalar_type>::matrix matrix;
	typedef math::state_functor<scalar_type> state_functor;

	/** Construct a scalar product functor from a given state functor and
	 * a vector of states that will be used for computing the scalar products.
	 */
	constant_state_matrix_functor(
			const matrix & A ) : my_A(A) {
	}
	;

	/** Virtual desctructor for possible derived classes */
	virtual ~constant_state_matrix_functor() {
	}
	;

	/**
	 * returns my_A
	 */
	virtual matrix map(const state& x) const {
		matrix res(my_A);
		return res;
	}
	;
private:
	const matrix & my_A;
};


}

#endif /* STATE_MATRIX_FUNCTOR_UTILITY_H_ */
