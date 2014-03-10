/*
 * vdom_matrix_operators.hpp
 *
 *  Created on: Oct 27, 2010
 *      Author: frehse
 */

#ifndef VDOM_MATRIX_OPERATORS_HPP_
#define VDOM_MATRIX_OPERATORS_HPP_

#include "vdom_matrix.h"
#include "vdom_matrix_utility.h"
#include "vdom_vector.h"

namespace math {

template<typename scalar_type> bool operator==(
		const vdom_matrix<scalar_type>& m1, const vdom_matrix<scalar_type>& m2) {
	if (m1.domain() != m2.domain() || m1.codomain() != m2.codomain()) {
		return false;
	}
	for (unsigned int i = 0; i < m1.size1(); i++) {
		for (unsigned int j = 0; j < m1.size2(); j++) {
			if (m1(i, j) != m2(i, j)) {
				return false;
			}
		}
	}
	return true;
}

template<typename scalar_type> bool operator!=(
		const vdom_matrix<scalar_type>& m1, const vdom_matrix<scalar_type>& m2) {
	return !(m1 == m2);
}

template<typename scalar_type> vdom_matrix<scalar_type> operator-(
		const vdom_matrix<scalar_type>& M) {
	return vdom_matrix<scalar_type> (M.codomain(), M.domain(),
			-M.get_matrix());
}


template<typename scalar_type> vdom_matrix<scalar_type> operator*(
		const vdom_matrix<scalar_type>& M1, const vdom_matrix<scalar_type>& M2) {
	if (M1.domain() == M2.codomain()) {
		return vdom_matrix<scalar_type> (M1.codomain(), M2.domain(),
				M1.get_matrix() * M2.get_matrix());
	} else {
		throw std::runtime_error(
				"missing implementation of matrix mult with remapping");
		return vdom_matrix<scalar_type> ();
	}
}

template<typename scalar_type> vdom_matrix<scalar_type> operator+(
		const vdom_matrix<scalar_type>& M1, const vdom_matrix<scalar_type>& M2) {
	if (M1.domain() == M2.codomain()) {
		return vdom_matrix<scalar_type> (M1.codomain(), M2.domain(),
				M1.get_matrix() + M2.get_matrix());
	} else {
		throw std::runtime_error(
				"missing implementation of matrix addition with remapping");
		return vdom_matrix<scalar_type> ();
	}
}

template<typename scalar_type> vdom_matrix<scalar_type> operator-(
		const vdom_matrix<scalar_type>& M1, const vdom_matrix<scalar_type>& M2) {
	if (M1.domain() == M2.codomain()) {
		return vdom_matrix<scalar_type> (M1.codomain(), M2.domain(),
				M1.get_matrix() - M2.get_matrix());
	} else {
		throw std::runtime_error(
				"missing implementation of matrix subtraction with remapping");
		return vdom_matrix<scalar_type> ();
	}
}

template<typename scalar_type> vdom_matrix<scalar_type> operator*(
		const scalar_type& s, const vdom_matrix<scalar_type>& M) {
	return vdom_matrix<scalar_type> (M.codomain(), M.domain(), s
			* M.get_matrix());
}

template<typename scalar_type> vdom_vector<scalar_type> operator*(
		const vdom_matrix<scalar_type>& M, const vdom_vector<scalar_type>& v) {
	return M.multiply_vector(v);
}

template<typename scalar_type> vdom_vector<scalar_type> operator*(
		const vdom_vector<scalar_type>& v, const vdom_matrix<scalar_type>& M) {
	if (M.codomain() == v.domain()) {
		return vdom_vector<scalar_type> (M.domain(),
				v.get_vector() * M.get_matrix());
	} else {
		return M.transpose() * v;
	}
}

template<typename scalar_type> vdom_matrix<scalar_type> matrix_exponential(
		const vdom_matrix<scalar_type>& A) {
	vdom_matrix<scalar_type> res(A.codomain(), A.domain(), matrix_exponential(
			A.get_matrix()));
	return res;
}

}

#endif /* VDOM_MATRIX_OPERATORS_HPP_ */
