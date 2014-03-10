#ifndef MATRIX_OPERATORS_H_
#define MATRIX_OPERATORS_H_

/**
 * @file
 * Mathematical operators for the matrix class.
 *
 * Mathematical operators have a particular mathematical
 * significance (unlike utility functions).
 */

#include "math/vector.h"
#include "math/matrix.h"

namespace math {

/** Element-wise matrix comparison.
 *
 * Returns true if two given matrices are element-wise equal.
 * Uses the standard != operator on elements (not a numerically approximative
 * comparison).
 *
 * If the dimensions of the matrices differ, the matrices
 * are considered not equal.
 */
template<typename scalar_type> bool operator==(const matrix<scalar_type>& m1,
		const matrix<scalar_type>& m2);

/** Element-wise matrix negative comparison.
 *
 * Returns true if two given matrices are not element-wise equal.
 * Uses the standard != operator on elements (not a numerically approximative
 * comparison).
 *
 * If the dimensions of the matrices differ, the matrices
 * are considered not equal.
 */
template<typename scalar_type> bool operator!=(const matrix<scalar_type>& m1,
		const matrix<scalar_type>& m2);

/** Element-wise sum of matrices with identical dimensions. */
template<typename scalar_type> matrix<scalar_type> operator+(
		const matrix<scalar_type>& M1, const matrix<scalar_type>& M2);

/** Element-wise difference of matrices with identical dimensions. */
template<typename scalar_type> matrix<scalar_type> operator-(
		const matrix<scalar_type>& M1, const matrix<scalar_type>& M2);

/** Matrix negation */
template<typename scalar_type> matrix<scalar_type> operator-(
		const matrix<scalar_type>& M);

/** Element-wise product of a matrix with a scalar. */
template<typename scalar_type> matrix<scalar_type> operator*(
		const scalar_type& s, const matrix<scalar_type>& M);

/** Element-wise product of a scalar with a matrix. */
template<typename scalar_type> matrix<scalar_type> operator*(
		const matrix<scalar_type>& M, const scalar_type& s);

/** Matrix-vector product. */
template<typename scalar_type> vector<scalar_type> operator*(const matrix<
		scalar_type>& M, const vector<scalar_type>& v);

/** Vector-matrix product. */
template<typename scalar_type> vector<scalar_type> operator*(const vector<
		scalar_type>& v, const matrix<scalar_type>& M);

/** Matrix-matrix product. */
template<typename scalar_type> matrix<scalar_type> operator*(const matrix<
		scalar_type>& M1, const matrix<scalar_type>& M2);

/** Matrix determinant.
 *
 * @author Matthias Althoff, 07 may 2010
 * rewritten since it should return a scalar */
template<typename scalar_type> scalar_type matrix_determinant(
		const matrix<scalar_type>& A);

/**
 * Matrix exponential, i.e., e^A.
 */
template<typename scalar_type> matrix<scalar_type> matrix_exponential(
		const matrix<scalar_type>& A);

/**
 * Matrix exponential, i.e., e^A, for Rational matrices.
 *
 * Specialization to Rational data type, which converts to double for the
 * exponentiation, and converts the result back to Rational.
 */
template<> inline matrix<Rational> matrix_exponential(const matrix<Rational>& A);

/** Compute the element-wise absolute of a matrix. */
template<typename scalar_type> matrix<scalar_type> absolute(const matrix<scalar_type>& A);

/** Compute a set of three special matrices used in ODE solving.
 *
 * They are given by the following equations:
 * - \f$ M1 = e^{A\delta} \f$
 * - \f$ M2 = \sum_{i=0}^{\infty} \frac{\delta^{i+1}}{(i+1)!} A^i \f$
 * - \f$ M3 = \sum_{i=0}^{\infty} \frac{\delta^{i+2}}{(i+2)!} A^i \f$
 *
 * */
template<typename scalar_type> void get_special_matrices(const matrix<scalar_type>& A, const scalar_type& delta,
		matrix<scalar_type>& special_matrix1,
		matrix<scalar_type>& special_matrix2,
		matrix<scalar_type>& special_matrix3);

/**
 * Creates an augmented matrix [A|B] from the passed matrices A and B. The columns of matrix B are
 * augmented with matrix A
 * @param A Matrix
 * @param B Matrix to augment with A.
 * @return Augmented matrix [A|B]
 */
template<typename scalar_type> matrix<scalar_type> augment_cols(const matrix<scalar_type>& A, const matrix<scalar_type>& B);
/**
 * Creates an augmented matrix[A/B], i.e., the rows of matrix B are augmented with the rows of matrix A.
 *
 * @param A
 * @param B
 * @return Augmented matrix [A/B]
 */
template<typename scalar_type> matrix<scalar_type> augment_rows(const matrix<scalar_type>&A, const matrix<scalar_type>&B);
/**
 * Computes the eigenvalues of the passed matrix A
 *
 * @return D_real Vector of real parts of the eigenvalues
 * @return D_imag Vector of imaginary parts of the eigenvalues
 * @return V Matrix with the eigenvectors as columns.
 */
template<typename scalar_type> void compute_eigenvalues(const matrix<scalar_type>& A,
																	   vector<scalar_type>& D_real,
																	   vector<scalar_type>& D_imag,
																	   matrix<scalar_type>& V);
/**
 * Computes the eigenvalues of the passed matrix A with Rational element types
 *
 * @param A Matrix of which the eigenvalues has to be computed.
 * @return D_real Vector of real parts of the eigenvalues
 * @return D_imag Vector of imaginary parts of the eigenvalues
 * @return V Matrix with the eigenvectors as columns.
 */
template<> inline void compute_eigenvalues(const matrix<Rational>&A,
										   vector<Rational>& D_real,
										   vector<Rational>& D_imag,
										   matrix<Rational>& V);

/**
 *
 * @param A
 * @param M
 * @param M_inv
 * @param A_min
 * @param eps
 */
template<typename scalar_type> void minimize_norm(const matrix<scalar_type>&A,
		matrix<scalar_type>&M, matrix<scalar_type>&M_inv,
		matrix<scalar_type>& A_min,
		const scalar_type eps);

/**
 *
 * @param A
 * @param M
 * @param M_inv
 * @param A_min
 * @param eps
 */
template<> void minimize_norm(const matrix<Rational>&A,
		matrix<Rational>&M, matrix<Rational>&M_inv,
		matrix<Rational>& A_min,
		const Rational eps);

/** Brings M to redcued row echelon form
 *
 * Adapted from the pseudocode at http://rosettacode.org/wiki/Reduced_row_echelon_form.
 */
template<typename scalar_type> void reduced_row_echelon_form(matrix<scalar_type>& M);

}



#include "matrix_operators.hpp"

#endif /*MATRIX_OPERATORS_H_*/
