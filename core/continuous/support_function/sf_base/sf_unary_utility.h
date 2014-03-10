/*
 * sf_unary_utility.h
 *
 *  Created on: Oct 28, 2010
 *      Author: frehse
 */

#ifndef SF_UNARY_UTILITY_H_
#define SF_UNARY_UTILITY_H_

#include "sf_unary.h"
#include "sf_unary_template.h"

namespace continuous {
namespace support_function {

/** Returns the set U centered around the origin and the original center in center_vec.
 *
 *  The center is computed as the center of the bounding box.
 *  If U is empty, center_vec is zero and a null pointer is returned.
 *  The domain of center_vec is used to construct the positional_vdomains.
 */
template<typename scalar_type>
continuous::support_function_provider::ptr construct_centered_set(
		const continuous::support_function_provider::const_ptr& U,
		math::vdom_vector<scalar_type>& center_vec) {
	using namespace continuous;

	// Use the hyperbox class to get the center
	support_function_provider::ptr res = support_function_provider::ptr();
	if (!U || math::definitely(U->is_empty())) {
		center_vec = math::vdom_vector<scalar_type>(center_vec.domain(),
				scalar_type(0));
		// let res be the null pointer
	} else {
		hyperbox<scalar_type> U_box;
		U_box = compute_bounding_box<scalar_type> (*U);

		if (U_box.is_finite()) {
			typename hyperbox<scalar_type>::vdom_vector_type U_center_tmp =
					U_box.compute_finite_center();
			U_center_tmp.remap(center_vec.domain());
			center_vec = U_center_tmp;

			// translate U by -U_center
			math::affine_map<scalar_type> translate_map(-center_vec);

			res = support_function_provider::ptr(
					new support_function::sf_unary<scalar_type>(U,
							translate_map));
		} else {
			// Don't center, return a clone of U
			center_vec = math::vdom_vector<scalar_type>(center_vec.domain(),
					scalar_type(0));
			res = continuous::support_function_provider::ptr(U->clone());
		}
	}
	return res;
}
;

/** An easy to use interface for constructing the convex hull of a class that does not
 * subclass support_function_provider
 *
 * This can be used, e.g., for polyhedron_collection.
 * */
template<typename scalar_type, class implementor>
continuous::support_function_provider::ptr construct_convex_hull(
		const boost::shared_ptr<implementor>& U_ptr) {
	typename sf_unary_template<scalar_type,implementor>::ptr chull_ptr(new sf_unary_template<scalar_type,implementor>(U_ptr));
	return chull_ptr;
}

}
}

#endif /* SF_UNARY_UTILITY_H_ */
