#ifndef MATRIX_UTILITY_H_
#define MATRIX_UTILITY_H_

/**
 * @file
 * Utility functions for the matrix class.
 *
 * Utility functions allow one construct or modify
 * matrices without being of particular mathematical
 * significance (like operators).
 *
 */

#include "math/scalar_types/rational.h"
#include "math/ublas_utility/ublas_utility.h"
#include "math/ublas_utility/expm.h"
#include "math/basic_functions.h"
#include "math/vdom/index_to_index_bimap.h"
#include "math/matrix.h"

namespace math {

/** Creates a diagonal matrix with dimensions n1 and n2 where the elements on
 * the diagonal a(i,i) are d_val and all other are o_val. */
template<typename scalar_type> matrix<scalar_type> diagonal_matrix(
		typename matrix<scalar_type>::size_type n1,
		typename matrix<scalar_type>::size_type n2, const scalar_type& d_val,
		const scalar_type& o_val);

/** Creates a square diagonal matrix with dimensions n where the elements on the
 * diagonal are d_val and all other are scalar_type(0). */
template<typename scalar_type> matrix<scalar_type> diagonal_matrix(
		typename matrix<scalar_type>::size_type n, const scalar_type& d_val);

/** Diagonal matrix from vector
 *
 * Creates a diagonal matrix with dimensions n1 and n2 where the elements on
 * the diagonal a(i,i)=v[i]. */
template<typename scalar_type> matrix<scalar_type> diagonal_matrix(
		typename matrix<scalar_type>::size_type n1,
		typename matrix<scalar_type>::size_type n2,
		const vector<scalar_type>& v);

/** Returns true if the elements on the
 * diagonal are d_val and all other are scalar_type(0). */
template<typename scalar_type> bool is_diagonal_matrix(
		const matrix<scalar_type>& M, const scalar_type& d_val);

/** Remap the matrix A according to the map f.
 *
 * The result is a m by n matrix R such that R(f(i),f(j))=A(i,j).
 * If m or n are omitted or zero, the respective dimension of R is determined by the max of
 * f(i) of the respective indices.
 *
 * Rows and columns that are not in f(i) are not R.
 * */
template<typename scalar_type> matrix<scalar_type> remap(const matrix<
		scalar_type>& A, const index_to_index_bimap& f, typename matrix<scalar_type>::size_type m=0, typename matrix<scalar_type>::size_type n=0);

/** Remap the rows of a matrix.
 *
 * The new matrix is the size of the highest index that f maps to.
 * Any rows that are not in the domain of f are also not in
 * the result. */
template<typename scalar_type> matrix<scalar_type> map_rows(const matrix<
		scalar_type>& A, const index_to_index_bimap& f);

/** Bring almost zero values of M to zero */
template<typename scalar_type> void snap_to_zero(matrix<scalar_type>& M);

typedef matrix<Rational> rational_matrix;
typedef matrix<double> double_matrix;

}

#include "matrix_utility.hpp"

#endif /*MATRIX_UTILITY_H_*/
