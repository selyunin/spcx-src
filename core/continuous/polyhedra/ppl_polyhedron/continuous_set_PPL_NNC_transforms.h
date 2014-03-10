/*
 * continuous_set_PPL_NNC_transforms.h
 *
 *  Created on: Mar 30, 2010
 *      Author: frehse
 */

#ifndef CONTINUOUS_SET_PPL_NNC_TRANSFORMS_H_
#define CONTINUOUS_SET_PPL_NNC_TRANSFORMS_H_

#include "continuous_set_PPL_NNC.h"
#include "core/continuous/continuous_set_operator_implementations/compute_transformation.h"
//#include "core/continuous/continuous_set_transforms/reset_affine_transform.h"

namespace continuous {

template<> class transformer<ppl_polyhedron::continuous_set_PPL_NNC> {
public:
	typedef ppl_polyhedron::continuous_set_PPL_NNC T;
	typedef ppl_polyhedron::continuous_set_PPL_NNC::ptr ptr;
	static continuous_set::ptr compute_or_assign(T& c,
			const constant_bound_time_elapse_transform& t) {
		c.assign_transformation(t);
		return c.get_ptr();
	}
	;
	static continuous_set::ptr compute(const T& c,
			const constant_bound_time_elapse_transform& t) {
		ptr res = ptr(c.clone());
		return compute_or_assign(*res, t);
	}
	;
	static continuous_set::ptr compute_or_assign(T& c,
			const reset_affine_transform<global_types::rational_type>& t) {
		throw std::runtime_error(
				"missing implementation of reset_affine_transform");
		return continuous_set::ptr();
	}
	;
	static continuous_set::ptr compute(const T& c,
			const reset_affine_transform<global_types::rational_type>& t) {
		ptr res = ptr(c.clone());
		return compute_or_assign(*res, t);
	}
	;
	static continuous_set::ptr compute_or_assign(T& c,
			const reset_affine_transform<global_types::float_type>& t) {
		math::affine_map<global_types::rational_type> M(t.convert_to<
				global_types::rational_type> ());
		reset_affine_transform<global_types::rational_type> trat(M);
		return compute_or_assign(c, trat);
	}
	;
	static continuous_set::ptr compute(const T& c,
			const reset_affine_transform<global_types::float_type>& t) {
		ptr res = ptr(c.clone());
		return compute_or_assign(*res, t);
	}
	;
	static continuous_set::ptr compute_or_assign(T& c,
			const reset_function_transform& t) {
		c.assign_transformation(t);
		return c.get_ptr();
	}
	;
	static continuous_set::ptr compute(const T& c,
			const reset_function_transform& t) {
		ptr res = ptr(c.clone());
		return compute_or_assign(*res, t);
	}
};

}

#endif /* CONTINUOUS_SET_PPL_NNC_TRANSFORMS_H_ */
