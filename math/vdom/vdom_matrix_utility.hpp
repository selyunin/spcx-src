/*
 * vdom_matrix_utility.hpp
 *
 *  Created on: Oct 27, 2010
 *      Author: frehse
 */

#ifndef VDOM_MATRIX_UTILITY_HPP_
#define VDOM_MATRIX_UTILITY_HPP_

#include "vdom_matrix_utility.h"

#include "vdom_matrix.h"
#include "math/ublas_utility/ublas_utility.h"
#include "math/ublas_utility/expm.h"
#include "math/basic_functions.h"
#include "math/vdom/index_to_index_bimap.h"
#include "math/matrix_utility.h"

namespace math {

template<typename scalar_type> vdom_matrix<scalar_type> diagonal_matrix(
		const positional_vdomain& codom, const positional_vdomain& dom,
		const scalar_type& d_val, const scalar_type& o_val) {
	vdom_matrix<scalar_type> m(codom, dom, diagonal_matrix(codom.size(),
			dom.size(), d_val, o_val));
	return m;
}

template<typename scalar_type> vdom_matrix<scalar_type> diagonal_matrix(
		const positional_vdomain& dom, const scalar_type& d_val) {
	return diagonal_matrix(dom, dom, d_val, scalar_type(0));
}

template<typename scalar_type> vdom_matrix<scalar_type> diagonal_matrix(
		const positional_vdomain& codom, const positional_vdomain& dom,
		const vdom_vector<scalar_type>& v) {
	vdom_matrix<scalar_type> m(codom, dom, diagonal_matrix(codom.size(),
			dom.size(), v.get_vector()));
	return m;
}

template<typename scalar_type> bool is_diagonal_matrix(
		const vdom_matrix<scalar_type>& M, const scalar_type& d_val) {
	return is_diagonal_matrix(M.get_matrix(),d_val);
}

template<typename scalar_type, template<typename > class resolver> vdom_matrix<
		scalar_type> compose(const vdom_matrix<scalar_type>& A1,
		const vdom_matrix<scalar_type>& A2, position_map& f1, position_map& f2,
		position_map& g1, position_map& g2, const resolver<scalar_type>& resol) {
	positional_vdomain C, D; // codomain and domain of A

	C = compose(A1.codomain(), A2.codomain(), f1, f2);
	D = compose(A1.domain(), A2.domain(), g1, g2);

	matrix<scalar_type> A(C.size(), D.size());

	/* Copy elements of A1 */
	typedef typename vdom_matrix<scalar_type>::size_type size_type;
	for (size_type i = 0; i < A1.size1(); ++i) {
		size_type f_i = f1.get_map(i);
		for (size_type j = 0; j < A1.size2(); ++j) {
			size_type g_j = g1.get_map(j);
			A(f_i, g_j) = A1(i, j);
		}
	}

	/* Copy elements of A2 */
	for (size_type i = 0; i < A2.size1(); ++i) {
		size_type f_i = f2.get_map(i);
		// we need to test for conflicts if this row received values from A1
		bool test_row_conflict = f1.in_codomain(f_i);
		for (size_type j = 0; j < A2.size2(); ++j) {
			size_type g_j = g2.get_map(j);
			if (!test_row_conflict)
				A(f_i, g_j) = A2(i, j);
			else {
				A(f_i, g_j) = resol.resolve(A(f_i, g_j), A2(i,
						j));
			}
		}
	}

	/* Return the result */
	return vdom_matrix<scalar_type> (C, D, A);
}

template<typename scalar_type> vdom_matrix<scalar_type> compose(
		const vdom_matrix<scalar_type>& A1, const vdom_matrix<scalar_type>& A2,
		position_map& f1, position_map& f2, position_map& g1, position_map& g2) {
	return compose<scalar_type, throwing_resolver> (A1, A2, f1, f2, g1, g2);
}

template<typename scalar_type> vdom_matrix<scalar_type> compose(
		const vdom_matrix<scalar_type>& A1, const vdom_matrix<scalar_type>& A2) {
	position_map f1,f2,g1,g2;
	return compose<scalar_type, throwing_resolver> (A1, A2, f1, f2, g1, g2);
}

template<typename scalar_type> void decompose_with_codom(const vdom_matrix<
		scalar_type>& A, vdom_matrix<scalar_type>& A1,
		vdom_matrix<scalar_type>& A2) {
	decompose_on_dom(A.codomain(),A,A1,A2);
}

template<typename scalar_type> void decompose_on_domain(const positional_vdomain& A1_dom, const vdom_matrix<
		scalar_type>& A, vdom_matrix<scalar_type>& A1,
		vdom_matrix<scalar_type>& A2) {
	// get the variables of the domain of A that are not in dom.
	variable_id_set domvars = A.domain().get_variable_ids();
	variable_id_set A1domvars = A1_dom.get_variable_ids();

	variable_id_set A2domvars = domvars;
		set_difference_assign(A2domvars, A1domvars);

	const positional_vdomain& Acodom = A.codomain();
	const positional_vdomain& Adom = A.domain();
	positional_vdomain A2dom = positional_vdomain(A2domvars);

	// For every variable in the codomain of A,
	// assign the corresponding column of A1 respectively A2.
	A1 = vdom_matrix<scalar_type> (Acodom, A1_dom);
	A2 = vdom_matrix<scalar_type> (Acodom, A2dom);

	const unsigned int n = Adom.size();
	const unsigned int m = Acodom.size();

	for (unsigned int i = 0; i < n; ++i) {
		unsigned int column_index;
		variable column_var = Adom.get_variable(i);
		if (A1_dom.in_domain(column_var, column_index)) {
			// assign the column to A1
			// no reordering necessary
			// since the codom of A and A1 are identical
			for (unsigned int j = 0; j < m; ++j) {
				A1(j, column_index) = A(j, i);
			}
		} else {
			// assign the column to A2
			// no reordering necessary
			// since the codom of A and A2 are identical
			unsigned int A2column_index = A2dom.pos(column_var);
			for (unsigned int j = 0; j < m; ++j) {
				A2(j, A2column_index) = A(j, i);
			}
		}
	}
}

}

#endif /* VDOM_MATRIX_UTILITY_HPP_ */
