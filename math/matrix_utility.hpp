/*
 * matrix_utility.hpp
 *
 *  Created on: Oct 27, 2010
 *      Author: frehse
 */

#ifndef MATRIX_UTILITY_HPP_
#define MATRIX_UTILITY_HPP_

#include "matrix_utility.h"

#include "math/scalar_types/rational.h"
#include "math/ublas_utility/ublas_utility.h"
#include "math/ublas_utility/expm.h"
#include "math/basic_functions.h"
#include "math/vdom/index_to_index_bimap.h"
#include "math/matrix.h"

namespace math {

template<typename scalar_type> matrix<scalar_type> diagonal_matrix(
		typename matrix<scalar_type>::size_type n1,
		typename matrix<scalar_type>::size_type n2, const scalar_type& d_val,
		const scalar_type& o_val) {
	matrix<scalar_type> m(n1, n2);
	for (unsigned int j = 0; j < n2; j++) {
		for (unsigned int i = 0; i < n1; i++) {
			if (i == j)
				m(i, j) = d_val;
			else
				m(i, j) = o_val;
		}
	}
	return m;
}

template<typename scalar_type> bool is_diagonal_matrix(
		const matrix<scalar_type>& M, const scalar_type& d_val) {
	for (unsigned int j = 0; j < M.size2(); j++) {
		for (unsigned int i = 0; i < M.size1(); i++) {
			if (i == j && M(i, j) != d_val) {
				return false;
			} else if (M(i, j) != scalar_type(0)) {
				return false;
			}
		}
	}
	return true;
}

template<typename scalar_type> matrix<scalar_type> diagonal_matrix(
		typename matrix<scalar_type>::size_type n, const scalar_type& d_val) {
	return diagonal_matrix(n, n, d_val, scalar_type(0));
}

template<typename scalar_type> matrix<scalar_type> diagonal_matrix(
		typename matrix<scalar_type>::size_type n1,
		typename matrix<scalar_type>::size_type n2,
		const vector<scalar_type>& v) {
	matrix<scalar_type> m(n1, n2, scalar_type(0));
	for (unsigned int j = 0; j < std::min(n1, n2); j++) {
		m(j, j) = v[j];
	}
	return m;
}

template<typename scalar_type> matrix<scalar_type> remap(const matrix<
		scalar_type>& A, const index_to_index_bimap& f, typename matrix<scalar_type>::size_type m=0, typename matrix<scalar_type>::size_type n=0) {

	typename matrix<scalar_type>::size_type n1 = A.size1();
	typename matrix<scalar_type>::size_type n2 = A.size2();
	// Find the max row index
	typename matrix<scalar_type>::size_type m1 = m;
	typename matrix<scalar_type>::size_type m2 = n;
	if (m==0)
	for (typename matrix<scalar_type>::size_type i = 0; i < n1; i++) {
		if (f.in_domain(i) && f.get_map(i) > m1)
			m1 = f.get_map(i);
	}
	if (n==0)
	for (typename matrix<scalar_type>::size_type i = 0; i < n2; i++) {
		if (f.in_domain(i) && f.get_map(i) > m2)
			m2 = f.get_map(i);
	}
	matrix<scalar_type> R(m1, m2);
	typename matrix<scalar_type>::size_type i_R;
	typename matrix<scalar_type>::size_type j_R;
	for (typename matrix<scalar_type>::size_type j = 0; j < n2; j++) {
		if (f.in_domain(j)) {
			j_R = f.get_map(j);
			for (typename matrix<scalar_type>::size_type i = 0; i < n1; i++) {
				if (f.in_domain(i)) {
					i_R = f.get_map(i);
					R(i_R, j_R) = A(i, j);
				}
			}
		}
	}
	return R;
}

template<typename scalar_type> matrix<scalar_type> map_rows(const matrix<
		scalar_type>& A, const index_to_index_bimap& f) {
	matrix<scalar_type> M(f.max_in_codomain() + 1, A.size2());

	for (unsigned int i = 0; i < A.size1(); i++) {
		M.row_proxy(f.get_map(i)) = A.row_proxy(i);
	}

	return M;
}

template<typename scalar_type> void snap_to_zero(matrix<scalar_type>& M) {
	using namespace math::numeric;
	for (size_t i = 0; i < M.size1(); ++i) {
		for (size_t j = 0; j < M.size2(); ++j) {
			if (is_MEQ(M(i, j), scalar_type(0))) {
				M(i, j) = scalar_type(0);
			}
		}
	}
}

}

#endif /* MATRIX_UTILITY_HPP_ */
