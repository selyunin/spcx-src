/*
 * affine_map_operators.h
 *
 *  Created on: Nov 9, 2009
 *      Author: frehse
 */

#ifndef AFFINE_MAP_OPERATORS_H_
#define AFFINE_MAP_OPERATORS_H_

/**
 * @file
 * Mathematical operators for the affine_map class.
 *
 * Mathematical operators have a particular mathematical
 * significance (unlike utility functions).
 */

#include "math/vdom/affine_map.h"

namespace math {

/** Returns true if two given maps are element-wise equal.
 *
 * Uses the standard != operator on elements (not a numerically approximative
 * comparison).
 *
 * If the domains or codomains of the maps differ, the matrices
 * are considered not equal.
 */
template<typename scalar_type> bool operator==(
		const affine_map<scalar_type>& m1, const affine_map<scalar_type>& m2);

/** Returns true if two given maps are not element-wise equal.
 *
 * Uses the standard != operator on elements (not a numerically approximative
 * comparison).
 *
 * If the domains or codomains of the maps differ, the matrices
 * are considered not equal.
 */
template<typename scalar_type> bool operator!=(
		const affine_map<scalar_type>& m1, const affine_map<scalar_type>& m2);

/** Returns the concatenation of affine maps M1 and M2.
 *
 * If M1:z=A1*y+b1 and M2:y=A2*x+b2, the result
 * is the map applying first M2 and then M1:
 *    M:z=A1*A2*x+A1*b2+b1.
 */
template<typename scalar_type> affine_map<scalar_type> concatenate(
		const affine_map<scalar_type>& M1, const affine_map<scalar_type>& M2);

/** Returns the parallel composition of M1 and M2.
 *
 * The parallel composition of two affine maps M1(x)=A1*x+b1 and M2(y)=A2*y+b2,
 * where x and y are named vectors over possibly different sets of variables,
 * is the map M=A*z+b for which
 * M(v)=M1(v) for all variables v in the domain of x, and
 * M(w)=M2(w) for all variables w in the domain of y.
 *
 * If the two maps disagree on a coefficient, the conflict is handled
 * by the resolver class.
 */
template<typename scalar_type, template<typename > class resolver> affine_map<
		scalar_type> compose(const affine_map<scalar_type>& M1,
		const affine_map<scalar_type>& M2, const resolver<scalar_type>& resol = resolver<scalar_type>());

/** Returns the parallel composition of M1 and M2.
 *
 * The parallel composition of two affine maps M1(x)=A1*x+b1 and M2(y)=A2*y+b2,
 * where x and y are named vectors over possibly different sets of variables,
 * is the map M=A*z+b for which
 * M(v)=M1(v) for all variables v in the domain of x, and
 * M(w)=M2(w) for all variables w in the domain of y.
 *
 * If such a map doesn't exist, an exception is thrown.
 */
template<typename scalar_type> affine_map<scalar_type> compose(
		const affine_map<scalar_type>& M1, const affine_map<scalar_type>& M2);

/** Returns the map resulting form adding a const vector to the map */
template<typename scalar_type> affine_map<scalar_type> operator+(
		const affine_map<scalar_type>& M, const vdom_vector<scalar_type>& v);

}

#include "affine_map_operators.hpp"

#endif /* AFFINE_MAP_OPERATORS_H_ */
