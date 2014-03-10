/*
 * affine_map_operators.hpp
 *
 *  Created on: Oct 27, 2010
 *      Author: frehse
 */

#ifndef AFFINE_MAP_OPERATORS_HPP_
#define AFFINE_MAP_OPERATORS_HPP_

#include "math/vdom/affine_map_operators.h"

namespace math {

template<typename scalar_type> bool operator==(
		const affine_map<scalar_type>& m1, const affine_map<scalar_type>& m2) {
	return (m1.get_A()==m2.get_A()) && (m1.get_b()==m2.get_b());
	/*
	if (m1.domain() == m2.domain() && m1.codomain() == m2.codomain()) {
		return (m1.get_A().get_matrix() == m2.get_A().get_matrix())
				&& (m1.get_b().get_vector() == m2.get_b().get_vector());
	} else {
		// return false if the variables of m2 aren't contained in those of m1
		if (m1.domain().contains_variables(m2.domain())
				&& m1.codomain().contains_variables(m2.codomain())) {
			affine_map<scalar_type> m2mapped(m2);
			m2mapped.reorder(m1.codomain(), m1.domain());
			return (m1.get_A().get_matrix() == m2mapped.get_A().get_matrix())
					&& (m1.get_b().get_vector() == m2mapped.get_b().get_vector());
		} else {
			return false;
		}
	}
	*/
}


template<typename scalar_type> bool operator!=(
		const affine_map<scalar_type>& m1, const affine_map<scalar_type>& m2) {
	return !(m1 == m2);
}

template<typename scalar_type> affine_map<scalar_type> concatenate(
		const affine_map<scalar_type>& M1, const affine_map<scalar_type>& M2) {
	if (M1.is_translation_coded() && M2.is_translation_coded()) {
		vdom_vector<scalar_type> b = M1.get_b() + M2.get_b();
		affine_map<scalar_type> M = affine_map<scalar_type> (b.domain(),b);
		return M;
	}
	if (!M1.is_translation_coded() && M2.is_translation_coded() && M1.codomain()==M2.codomain()) {
		affine_map<scalar_type> M = affine_map<scalar_type> (
				M1.get_A(),M1.get_b() + M2.get_b());
		return M;
	}
	if (M1.is_translation_coded() && !M2.is_translation_coded() && M1.codomain()==M2.codomain()) {
		affine_map<scalar_type> M = affine_map<scalar_type> (
				M2.get_A(),M1.get_b() + M2.get_b());
		return M;
	}

	affine_map<scalar_type> M(M1.get_A() * M2.get_A(),
			M1.get_A() * M2.get_b() + M1.get_b());
	return M;

}

template<typename scalar_type, template<typename > class resolver> affine_map<
		scalar_type> compose(const affine_map<scalar_type>& M1,
		const affine_map<scalar_type>& M2, const resolver<scalar_type>& resol) {

//	std::cout << "composing " << M1 << " with " << M2 << std::endl << std::flush;

	// If the maps are void, return a void map
	if (M1.is_void() || M2.is_void()) {
		positional_vdomain cdom = compose(M1.domain(), M2.domain());
		return math::affine_map<scalar_type>::void_map(cdom);
	}

	position_map f1;
	position_map f2;
	affine_map<scalar_type> res;
	if (M1.is_translation_coded() && M2.is_translation_coded()) {
//		std::cout << "b1:" << M1.get_b() << ", b2:" << M2.get_b() << std::endl;

		positional_vdomain C, D; // codomain and domain of A
		C = compose(M1.codomain(), M2.codomain(), f1, f2);
		vector<scalar_type> b_vec;
		b_vec = compose<scalar_type, resolver> (M1.get_b().get_vector(),
				M2.get_b().get_vector(), f1, f2, resol);
		vdom_vector<scalar_type> b(C, b_vec);

		res = affine_map<scalar_type> (C,b);
	} else {
		position_map g1;
		position_map g2;

		// Compose A1 and A2 into A
		vdom_matrix<scalar_type> A;
		A = compose<scalar_type, resolver> (M1.get_A(), M2.get_A(), f1, f2, g1,
				g2);

		// Compose b1 and b2 into b using the maps from A
		vector<scalar_type> b_vec;
		b_vec = compose<scalar_type, resolver> (M1.get_b().get_vector(),
				M2.get_b().get_vector(), f1, f2, resol);
		vdom_vector<scalar_type> b(A.codomain(), b_vec);

		res = affine_map<scalar_type> (A, b);
	}
//	std::cout << "result :" << res <<std::endl;
	return res;
}

template<typename scalar_type> affine_map<scalar_type> compose(
		const affine_map<scalar_type>& M1, const affine_map<scalar_type>& M2) {
	return compose<scalar_type, throwing_resolver> (M1, M2);
}

template<typename scalar_type> affine_map<scalar_type> operator+(
		const affine_map<scalar_type>& M, const vdom_vector<scalar_type>& v) {
	affine_map<scalar_type> K = M;
	K+=v;
	return K;
}


}

#endif /* AFFINE_MAP_OPERATORS_HPP_ */
