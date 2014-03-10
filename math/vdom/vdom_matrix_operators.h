#ifndef vdom_MATRIX_OPERATORS_H_
#define vdom_MATRIX_OPERATORS_H_

/**
 * @file
 * Mathematical operators for the vdom_matrix class.
 *
 * Mathematical operators have a particular mathematical
 * significance (unlike utility functions).
 */

#include "vdom_matrix.h"
#include "vdom_matrix_utility.h"
#include "vdom_vector.h"

namespace math {

/** Returns true if two given matrices are element-wise equal.
 *
 * Uses the standard != operator on elements (not a numerically approximative
 * comparison).
 *
 * If the domains or codomains of the matrices differ, the matrices
 * are considered not equal.
 */
template<typename scalar_type> bool operator==(
		const vdom_matrix<scalar_type>& m1, const vdom_matrix<scalar_type>& m2);

/** Returns true if two given matrices are not element-wise equal.
 *
 * Uses the standard != operator on elements (not a numerically approximative
 * comparison).
 *
 * If the domains or codomains of the matrices differ, the matrices
 * are considered not equal.
 */
template<typename scalar_type> bool operator!=(
		const vdom_matrix<scalar_type>& m1, const vdom_matrix<scalar_type>& m2);

/** Matrix negation */
template<typename scalar_type> vdom_matrix<scalar_type> operator-(
		const vdom_matrix<scalar_type>& M);

/** Matrix-matrix product.
 *
 * Throws if the domain of M1 is not equal to the codomain of M2.
 * @todo Treat the cases where this could be resolved by reordering. */
template<typename scalar_type> vdom_matrix<scalar_type> operator*(
		const vdom_matrix<scalar_type>& M1, const vdom_matrix<scalar_type>& M2);

/** Scalar-matrix product. */
template<typename scalar_type> vdom_matrix<scalar_type> operator*(
		const scalar_type& s, const vdom_matrix<scalar_type>& M);

/** Matrix-vector product. */
template<typename scalar_type> vdom_vector<scalar_type> operator*(
		const vdom_matrix<scalar_type>& M, const vdom_vector<scalar_type>& v);

/** Vector-matrix product. */
template<typename scalar_type> vdom_vector<scalar_type> operator*(
		const vdom_vector<scalar_type>& v, const vdom_matrix<scalar_type>& M);

/** Computes Euler exponentiation of a square matrix, i.e. e^A.
 */
template<typename scalar_type> vdom_matrix<scalar_type> matrix_exponential(
		const vdom_matrix<scalar_type>& A);

}

#include "vdom_matrix_operators.hpp"

#endif /*vdom_MATRIX_OPERATORS_H_*/
