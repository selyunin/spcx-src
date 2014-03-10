/*
 * continuous_set_collection_operators.h
 *
 *  Created on: Mar 30, 2010
 *      Author: frehse
 */

#ifndef CONTINUOUS_SET_COLLECTION_OPERATORS_H_
#define CONTINUOUS_SET_COLLECTION_OPERATORS_H_

#include "continuous_set_collection.h"
#include "continuous_set_operator_implementations/compute_transformation.h"

namespace continuous {

template<> class transformer<continuous_set_collection> {
public:
	static continuous_set::ptr compute(const continuous_set_collection& c,
			const constant_bound_time_elapse_transform& t);

	static continuous_set::ptr compute_or_assign(continuous_set_collection& c,
			const constant_bound_time_elapse_transform& t);

	static continuous_set::ptr compute(const continuous_set_collection& c,
			const reset_affine_transform<global_types::rational_type>& t);

	static continuous_set::ptr compute(const continuous_set_collection& c,
			const reset_affine_transform<global_types::float_type>& t);

	static continuous_set::ptr compute_or_assign(continuous_set_collection& c,
			const reset_affine_transform<global_types::rational_type>& t);

	static continuous_set::ptr compute_or_assign(continuous_set_collection& c,
			const reset_affine_transform<global_types::float_type>& t);

	static continuous_set::ptr compute(const continuous_set_collection& c,
			const reset_function_transform& t);

	static continuous_set::ptr compute_or_assign(continuous_set_collection& c,
			const reset_function_transform& t);
};

}

#endif /* CONTINUOUS_SET_COLLECTION_OPERATORS_H_ */
