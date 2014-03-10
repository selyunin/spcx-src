/*
 * constr_polyhedron_transforms.h
 *
 *  Created on: Mar 30, 2010
 *      Author: frehse
 */

#ifndef CONSTR_POLYHEDRON_TRANSFORMS_H_
#define CONSTR_POLYHEDRON_TRANSFORMS_H_

#include "constr_polyhedron_operators.h"
#include "core/continuous/continuous_set_operator_implementations/compute_transformation.h"
#include "core/continuous/continuous_set_transforms/reset_affine_transform.h"

namespace continuous {

template<typename scalar_type> class transformer<constr_polyhedron<scalar_type> > {
public:
	typedef constr_polyhedron<scalar_type> T;
	static continuous_set::ptr compute(const T& c,
			const constant_bound_time_elapse_transform& t) {
		throw std::runtime_error(
				"missing implementation of constant_bound_time_elapse_transform");
		return continuous_set::ptr();
	}
	;
	static continuous_set::ptr compute_or_assign(T& c,
			const constant_bound_time_elapse_transform& t) {
		return compute(c, t);
	}
	;
	static continuous_set::ptr compute(const T& c,
			const reset_affine_transform<global_types::rational_type>& t) {
		math::affine_map<scalar_type> M(t.template convert_to<
				scalar_type> ());
		return continuous_set::ptr(new constr_polyhedron<scalar_type> (
				apply_map(c, M)));
	}
	;

	static continuous_set::ptr compute_or_assign(T& c,
			const reset_affine_transform<global_types::rational_type>& t) {
		return compute(c, t);
	}
	;
	static continuous_set::ptr compute(const T& c,
			const reset_affine_transform<global_types::float_type>& t) {
		math::affine_map<scalar_type> M(t.template convert_to<
				scalar_type> ());
		return continuous_set::ptr(new constr_polyhedron<scalar_type> (
				apply_map(c, M)));
	}
	;

	static continuous_set::ptr compute_or_assign(T& c,
			const reset_affine_transform<global_types::float_type>& t) {
		return compute(c, t);
	}
	;

	static continuous_set::ptr compute(const T& c,
			const reset_function_transform& t) {
		throw std::runtime_error(
				"missing implementation of reset_function_transform");
		return continuous_set::ptr();
	}
	static continuous_set::ptr compute_or_assign(T& c,
			const reset_function_transform& t) {
		return compute(c, t);
	}
	;
};

}
#endif /* CONSTR_POLYHEDRON_TRANSFORMS_H_ */
