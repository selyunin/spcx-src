#ifndef vdom_MATRIX_UTILITY_H_
#define vdom_MATRIX_UTILITY_H_

/**
 * @file
 * Utility functions for the vdom_matrix class.
 *
 * Utility functions allow one construct or modify
 * matrices without being of particular mathematical
 * significance (like operators).
 *
 */

#include "vdom_matrix.h"
#include "math/matrix_utility.h"

namespace math {

/** Creates a diagonal matrix with dimensions n1 and n2 where the elements on
 * the diagonal a(i,i) are d_val and all other are o_val. */
template<typename scalar_type> vdom_matrix<scalar_type> diagonal_matrix(
		const positional_vdomain& codom, const positional_vdomain& dom,
		const scalar_type& d_val, const scalar_type& o_val);

/** Creates a square diagonal matrix with dimensions n where the elements on the
 * diagonal are d_val and all other are scalar_type(0). */
template<typename scalar_type> vdom_matrix<scalar_type> diagonal_matrix(
		const positional_vdomain& dom, const scalar_type& d_val);

/** Diagonal matrix from vector
 *
 * Creates a diagonal matrix with dimensions n1 and n2 where the elements on
 * the diagonal a(i,i)=v[i]. */
template<typename scalar_type> vdom_matrix<scalar_type> diagonal_matrix(
		const positional_vdomain& codom, const positional_vdomain& dom,
		const vdom_vector<scalar_type>& v);

/** Returns true if the elements on the
 * diagonal are d_val and all other are scalar_type(0). */
template<typename scalar_type> bool is_diagonal_matrix(
		const matrix<scalar_type>& M, const scalar_type& d_val);

/** Map two matrices to a common domain using a conflict resolver
 * and returns the resulting matrix A.
 *
 * The conflict resolver determines what happens if the matrices have
 * elements that get mapped to the same row in the result, and
 * these elements are different from each other.
 * The resolver is templated with scalar_type.
 *
 * f1 and f2 map the codomain of A1 and A2 to the resulting matrix.
 * g1 and g2 map the domain of A1 and A2 to the resulting matrix.
 *
 * If there is no conflict, then for all (i,j) in the index space of A1,
 *    A(f1(i),g1(j))=A1(i,j),
 * and similarly for A2.
 * Even stronger, for any y=A*x, the projection onto the variables of
 * A1 satisfies y1=A1*x1 and similarly y2=A2*x2.
 * */
template<typename scalar_type, template<typename > class resolver> vdom_matrix<
		scalar_type> compose(const vdom_matrix<scalar_type>& A1,
		const vdom_matrix<scalar_type>& A2, position_map& f1, position_map& f2,
		position_map& g1, position_map& g2, const resolver<scalar_type>& resol = resolver<scalar_type>());

/** Map two matrices to a common domain and returns the resulting matrix A.
 *
 * Throws if the matrices have elements that get mapped to the same row
 * in the result, and these elements are different from each other.
 *
 * f1 and f2 map the codomain of A1 and A2 to the resulting matrix.
 * g1 and g2 map the domain of A1 and A2 to the resulting matrix.
 *
 * If there is no conflict, then for all (i,j) in the index space of A1,
 *    A(f1(i),g1(j))=A1(i,j),
 * and similarly for A2.
 * Even stronger, for any y=A*x, the projection onto the variables of
 * A1 satisfies y1=A1*x1 and similarly y2=A2*x2.
 * */
template<typename scalar_type> vdom_matrix<scalar_type> compose(
		const vdom_matrix<scalar_type>& A1, const vdom_matrix<scalar_type>& A2,
		position_map& f1, position_map& f2, position_map& g1, position_map& g2);

/** Map two matrices to a common domain and returns the resulting matrix A.
 *
 * Throws if the matrices have elements that get mapped to the same row
 * in the result, and these elements are different from each other.
 *
 * f1 and f2 map the codomain of A1 and A2 to the resulting matrix.
 * g1 and g2 map the domain of A1 and A2 to the resulting matrix.
 *
 * If there is no conflict, then for all (i,j) in the index space of A1,
 *    A(f1(i),g1(j))=A1(i,j),
 * and similarly for A2.
 * Even stronger, for any y=A*x, the projection onto the variables of
 * A1 satisfies y1=A1*x1 and similarly y2=A2*x2.
 * */
template<typename scalar_type> vdom_matrix<scalar_type> compose(
		const vdom_matrix<scalar_type>& A1, const vdom_matrix<scalar_type>& A2);

/** Decompose a matrix into its codomain plus a domain consisting of
 * the remaining variables.
 *
 * Returns matrices A1 and A2 such that
 * - [A1 A2] is equivalent to A up to reordering variables, and
 * - A1.domain()==A1.codomain(),
 * - A2.codomain()==A1.codomain().
 *
 * If the domain of A does not contain all variables of the codomain
 * of A, the corresponding entries are filled with zeros.
 * This preserves the semantics of a map: [A1 A2]*x == A*x for any x.
 *
 * @remark A2 (and its domains) can be empty.
 */
template<typename scalar_type> void decompose_with_codom(const vdom_matrix<
		scalar_type>& A, vdom_matrix<scalar_type>& A1,
		vdom_matrix<scalar_type>& A2);

/** Decompose a matrix into a given domain and the remaining variables.
 *
 * Returns matrices A1 and A2 such that
 * - [A1 A2] is equivalent to A up to reordering variables, and
 * - A1.domain()==A1_dom,
 * - A2.codomain()==A1.codomain()==A.codomain().
 *
 * If the domain of A does not contain all variables of A1_dom
 * of A, the corresponding entries are filled with zeros.
 * This preserves the semantics of a map: [A1 A2]*x == A*x for any x.
 *
 * @remark A2 (and its domains) can be empty.
 */
template<typename scalar_type> void decompose_on_domain(
		const positional_vdomain& A1_dom, const vdom_matrix<scalar_type>& A,
		vdom_matrix<scalar_type>& A1, vdom_matrix<scalar_type>& A2);

}

#include "vdom_matrix_utility.hpp"

#endif /*vdom_MATRIX_UTILITY_H_*/
