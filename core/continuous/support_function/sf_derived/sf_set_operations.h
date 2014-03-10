/*
 * sf_set_operations.h
 *
 *  Created on: Nov 12, 2012
 *      Author: notroot
 */

#ifndef SF_SET_OPERATIONS_H_
#define SF_SET_OPERATIONS_H_

/** This file contains miscellaneous operations on sf_set objects */

#include "sf_sum.h"

namespace continuous {
namespace support_function {

/** Returns a pointer to the set sfp mapped with M and bloated by W
 *
 * Formally, if sfp represents the set X, then the returned set is
 * Y = MX + W,
 * where addition is Minkowski sum.
 */
template<typename scalar_type>
typename sf_set<scalar_type>::ptr map_and_bloat(
		const support_function_provider::const_ptr& sfp,
		const math::affine_map<scalar_type>& M,
		const support_function_provider::const_ptr& W) {

	// the result is the minkowski sum of the mapped states and the mapped inputs
	// get the mapped states
	typename sf_set<scalar_type>::ptr mapped_states(
			new support_function::sf_unary<scalar_type>(sfp, M));
	// form the minkowski sum with the inputs
	typename sf_set<scalar_type>::ptr mapped_sf;
	if (W) {
		mapped_sf = typename sf_set<scalar_type>::ptr(
				new sf_sum<scalar_type>(mapped_states, W));
	} else {
		mapped_sf = mapped_states;
	}
	return mapped_sf;
}
;

/** Returns a pointer to the set sfp  bloated by W
 *
 * Formally, if sfp represents the set X, then the returned set is
 * Y = X + W,
 * where addition is Minkowski sum.
 */
template<typename scalar_type>
support_function_provider::const_ptr bloat(
		const support_function_provider::const_ptr& sfp,
		const support_function_provider::const_ptr& W) {

	// form the minkowski sum with the inputs
	support_function_provider::const_ptr mapped_sf;
	if (W) {
		mapped_sf = typename sf_set<scalar_type>::ptr(
				new sf_sum<scalar_type>(sfp, W));
	} else {
		mapped_sf = sfp;
	}
	return mapped_sf;
}
;

}
}

#endif /* SF_SET_OPERATIONS_H_ */
