/*
 * affine_map_utility.hpp
 *
 *  Created on: Oct 27, 2010
 *      Author: frehse
 */

#ifndef AFFINE_MAP_UTILITY_HPP_
#define AFFINE_MAP_UTILITY_HPP_

#include "affine_map_utility.h"

#include "math/vdom/affine_map.h"
#include "math/vdom/affine_map_operators.h"

namespace math {

template<typename scalar_type>
void map_to_common_domain(vdom_matrix<scalar_type>& A, vdom_vector<
		scalar_type>& b) {
	if (A.codomain() != b.domain()) {
		position_map f1;
		position_map f2;
		positional_vdomain C;
		C = compose(A.codomain(), b.domain(), f1, f2);
		matrix<scalar_type> Am = map_rows(A.get_matrix(), f1);
		vector<scalar_type> bm = map(b.get_vector(), f2);
		A = vdom_matrix<scalar_type> (C, A.domain(), Am);
		b = vdom_vector<scalar_type> (C, bm);
	}
}

template<typename scalar_type>
void decompose_on_domain(const positional_vdomain& dom, const math::affine_map<scalar_type>& M,
		math::affine_map<scalar_type>& M1,
		math::affine_map<scalar_type>& M2) {
	if (M.domain() != dom) {
		// split A into square A1 and A2
		math::vdom_matrix<scalar_type> A1;
		math::vdom_matrix<scalar_type> A2;
		math::decompose_on_domain(dom, M.get_A(), A1, A2);

		M1=math::affine_map<scalar_type>(A1,M.get_b());
		M2=math::affine_map<scalar_type>(A2);
		LOGGER_OS(DEBUG7,__FUNCTION__) << "split map " << M << " into " << M1 << " and " << M2
							<< " with the latter domain " << M2.domain() << std::endl;
	} else {
		M1 = M;
		M2 = math::affine_map<scalar_type>();
	}
}

template<typename scalar_type>
void separate_states_from_inputs(const math::affine_map<scalar_type>& M,
		math::affine_map<scalar_type>& Mstate,
		math::affine_map<scalar_type>& Minput) {
	if (M.domain() != M.codomain()) {
		// split A into square A1 and A2
		decompose_on_domain(M.codomain(),M,Mstate,Minput);
		if (Mstate.domain()!=Mstate.codomain()) {
			Mstate.reorder(M.codomain(),M.codomain());
		}
//		LOGGER_OS(DEBUG7,__FUNCTION__) << "split map " << M << " into " << Mstate << " and " << Minput
//							<< " with input domain " << Minput.domain() << std::endl;
	} else {
		Mstate = M;
		Minput = math::affine_map<scalar_type>();
	}
}
}

#endif /* AFFINE_MAP_UTILITY_HPP_ */
