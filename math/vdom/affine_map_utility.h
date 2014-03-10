/*
 * affine_map_utility.h
 *
 *  Created on: Nov 9, 2009
 *      Author: frehse
 */

#ifndef AFFINE_MAP_UTILITY_H_
#define AFFINE_MAP_UTILITY_H_

/**
 * @file
 * Utility functions for the affine_map class.
 *
 * Utility functions allow one construct or modify
 * matrices without being of particular mathematical
 * significance (like operators).
 *
 */

#include "math/vdom/affine_map.h"
#include "math/vdom/affine_map_operators.h"

namespace math {

/** Map the codomain of A and the domain of b to a common domain. */
template<typename scalar_type>
void map_to_common_domain(vdom_matrix<scalar_type>& A,
		vdom_vector<scalar_type>& b);

/**
 * Separate state from inputs variables.
 *
 * Given ode affine dynamics
 *    x'=A0[x;u]+b0
 * separate into
 *    x'=Ax+b0+Bu.
 * with M=(A0,b0), Mstate=(A,b0), Minput=(B,0).
 * It holds
 * - Mstate.codomain()==Mstate.domain(),
 * - Mstate.codomain()==Minput.codomain(),
 * - Mstate.domain() and Minput.domain() have no common variables.
 */
template<typename scalar_type>
void separate_states_from_inputs(const math::affine_map<scalar_type>& M,
		math::affine_map<scalar_type>& Mstate,
		math::affine_map<scalar_type>& Minput);


/**
 * Separate M into a map over a given domain and a map over the rest of the variables
 *
 * Let the variables of dom be x, then given ode affine dynamics
 *    y'=A[x;u]+b
 * separate into
 *    y'=A1x+b+A2u.
 * with M=(A,b), M1=(A1,b), M2=(A2,0).
 * It holds
 * - M1.codomain()==M.codomain(),
 * - M1.domain()==dom,
 * - M1.domain() union M2.domain()==M.domain(),
 * - M1.domain() and M2.domain() have no common variables.
 */
template<typename scalar_type>
void decompose_on_domain(const positional_vdomain& dom, const math::affine_map<scalar_type>& M,
		math::affine_map<scalar_type>& M1,
		math::affine_map<scalar_type>& M2);

}

#include "affine_map_utility.hpp"

#endif /* AFFINE_MAP_UTILITY_H_ */
